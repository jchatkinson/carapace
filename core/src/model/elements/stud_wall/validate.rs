//! Layout and station validation, with the offending item named in the error.

use super::*;

pub(super) fn dedup_x(xs: &mut Vec<f64>, l: f64) {
    xs.sort_by(f64::total_cmp);
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-12 * l.max(1.));
}
pub(super) fn validate_layout<K: NodeKey>(
    l: &StudWallLayout,
    stations: &[WallStation<K>],
) -> Result<(), StudWallError> {
    for (name, v) in [
        ("L", l.length),
        ("H", l.height),
        ("top.ea", l.top_ea),
        ("bottom.ea", l.bottom_ea),
        ("bottom.ei", l.bottom_ei),
    ] {
        if !v.is_finite() || v <= 0. {
            return Err(StudWallError(format!("{name} must be positive and finite")));
        }
    }
    if let Some(z) = &l.floor {
        for (name, v, ok) in [
            ("floor.depth", z.depth, z.depth > 0.),
            ("floor.ga", z.ga, z.ga > 0.),
            ("floor.rim_ea", z.rim_ea, z.rim_ea >= 0.),
        ] {
            if !v.is_finite() || !ok {
                return Err(StudWallError(format!(
                    "{name} must be finite and {}",
                    if name == "floor.rim_ea" {
                        "nonnegative"
                    } else {
                        "positive"
                    }
                )));
            }
        }
    }
    if l.top_distribution == ChordDistribution::Beam
        && l.top_ei.is_none_or(|x| !x.is_finite() || x <= 0.)
    {
        return Err(StudWallError(
            "top.ei must be positive and finite for Beam".into(),
        ));
    }
    for (i, s) in l.studs.iter().enumerate() {
        if !s.x.is_finite() || !s.ea.is_finite() || s.ea <= 0. {
            return Err(StudWallError(format!(
                "studs[{i}] x and ea must be finite and ea positive"
            )));
        }
        if i > 0 && s.x <= l.studs[i - 1].x {
            return Err(StudWallError(format!(
                "studs[{i}] must be strictly increasing in x"
            )));
        }
    }
    if l.studs.first().is_none_or(|s| s.x != 0.)
        || l.studs
            .last()
            .is_none_or(|s| (s.x - l.length).abs() > 1e-10 * l.length)
    {
        return Err(StudWallError("studs must start at 0 and end at L".into()));
    }
    for (i, o) in l.openings.iter().enumerate() {
        if !(0. < o.x0
            && o.x0 < o.x1
            && o.x1 < l.length
            && 0. < o.head_height
            && o.head_height < l.height)
        {
            return Err(StudWallError(format!(
                "openings[{i}] has invalid bounds or head_height"
            )));
        }
        if !o.header_ei.is_finite()
            || o.header_ei <= 0.
            || !o.jack_ea.is_finite()
            || o.jack_ea <= 0.
        {
            return Err(StudWallError(format!(
                "openings[{i}] header_ei and jack_ea must be positive and finite"
            )));
        }
        // A stud row at `x0` or `x1` is a king stud: it shares the jack's bottom-plate node (one end post).
        if i > 0
            && (o.x0 <= l.openings[i - 1].x1
                || !l
                    .studs
                    .iter()
                    .any(|s| s.x >= l.openings[i - 1].x1 && s.x <= o.x0))
        {
            return Err(StudWallError(format!(
                "openings[{i}] overlaps or has no full-height stud separation"
            )));
        }
    }
    if let Some(s) = &l.sheathing {
        match &s.law {
            WallShearLaw::Linear { ga } if !ga.is_finite() || *ga <= 0. => {
                return Err(StudWallError(
                    "sheathing.ga must be positive and finite".into(),
                ));
            }
            WallShearLaw::O86 {
                bv,
                nail_spacing,
                nail_diameter,
            } if [bv, nail_spacing, nail_diameter]
                .iter()
                .any(|x| !x.is_finite() || **x <= 0.) =>
            {
                return Err(StudWallError(
                    "O86 bv, nail_spacing and nail_diameter must be positive and finite".into(),
                ));
            }
            WallShearLaw::Material(m)
                if !m.initial_tangent().is_finite() || m.initial_tangent() <= 0. =>
            {
                return Err(StudWallError(
                    "sheathing material must have a positive finite initial tangent".into(),
                ));
            }
            _ => {}
        }
    }
    if let Some(b) = &l.bearing {
        if l.foundation {
            return Err(StudWallError("bearing and foundation are exclusive".into()));
        }
        if !b.kv.is_finite() || b.kv <= 0. || !b.kh.is_finite() || b.kh <= 0. {
            return Err(StudWallError(
                "bearing.kv and bearing.kh must be positive and finite".into(),
            ));
        }
    }
    for (i, s) in stations.iter().enumerate() {
        if s.chord == WallChord::Base {
            if l.bearing.is_none() {
                return Err(StudWallError(format!(
                    "stations[{i}] is a base station without bearing"
                )));
            }
            if !s.x.is_finite() {
                return Err(StudWallError(format!("stations[{i}].x must be finite")));
            }
            continue;
        }
        if !s.x.is_finite() || s.x < 0. || s.x > l.length {
            return Err(StudWallError(format!("stations[{i}].x must be in [0,L]")));
        }
        if s.x == 0. || s.x == l.length {
            return Err(StudWallError(format!(
                "stations[{i}] duplicates a corner location"
            )));
        }
        if stations[..i]
            .iter()
            .any(|p| p.chord == s.chord && (p.x - s.x).abs() < 1e-12 * l.length)
        {
            return Err(StudWallError(format!(
                "stations[{i}] duplicates another station location"
            )));
        }
        if s.chord == WallChord::Bottom
            && l.openings
                .iter()
                .any(|o| o.door && s.x > o.x0 && s.x < o.x1)
        {
            return Err(StudWallError(format!("stations[{i}] lies inside a door")));
        }
        // A top station between studs (cripples count) needs the plate's EI: with `Direct` the plate spans between the
        // studs either side as a simply supported beam, so its load reaches them by the lever rule.
        if s.chord == WallChord::Top
            && l.top_distribution == ChordDistribution::Direct
            && !l
                .studs
                .iter()
                .any(|w| (w.x - s.x).abs() < 1e-12 * l.length.max(1.))
            && l.top_ei.is_none_or(|x| !x.is_finite() || x <= 0.)
        {
            return Err(StudWallError(format!(
                "stations[{i}] is a top station between studs: top.ei must be positive and finite (the plate spans between the studs)"
            )));
        }
    }
    Ok(())
}
