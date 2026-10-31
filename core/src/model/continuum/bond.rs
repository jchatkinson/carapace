//! Bond-slip discretization of an embedded bar (plan `calcs.app/docs/concrete/csfm/bond-plan.md` §4; CSFM book §3.3.3, §3.6.5).
//!
//! [`bond_polyline`] takes the host-bound segments `embed_polyline_in` produced and replaces the perfect-bond kinematics with
//! slipping ones. At every bar point `k`:
//!
//! - a **concrete image** node `C_k` with no DOFs of its own: `u_C = sum_i N_i(x_k) u_host,i` (two linear constraints, the host
//!   being the element of the segment that starts at the point, or ends there for the last point);
//! - a **rod node** `R_k` with free `ux, uy`. Its normal motion is tied to `C_k` by a linear constraint, so the one free DOF
//!   is the slip along the tangent;
//! - an oriented `ZeroLength` between them whose local x is the tangent: the bond spring, a `Steel01` force law
//!   (`e0 = Gb P l_k`, `fy = fbd P l_k`, hardening `Rf`) with `l_k` the tributary length. It reports `BondSlip` utilization;
//! - at the bar ends, an optional parallel anchorage spring (`Fau`, `Ka`, hardening `KSN`) reporting `AnchorageSlip`.
//!
//! The rod segments become [`EmbeddedBar`]s between rod nodes (the `Rod` host), so the bar law (`CsfmRebar`, rupture, tension
//! stiffening routing, recorders) is unchanged. A vertex whose turn angle exceeds `lock_angle` is **slip-locked**: no rod node,
//! the segments attach to `C_k` and no bond spring acts there (a bend as a bearing anchorage, book §3.6.5).
//!
//! Units follow the model (N, mm, MPa): `perimeter` mm, `fbd` MPa, `gb` MPa/mm, slips mm.

use crate::model::{
    Domain, EmbeddedBar, EmbeddedBarHost, LimitKind, Material, ModelError, Node, NodeId,
    Orientation, ZeroLength,
};

/// An anchorage end spring: yield force `fau` (N), elastic stiffness `ka` (N/mm), plastic slope `ksn * ka`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EndSpring {
    pub fau: f64,
    pub ka: f64,
    pub ksn: f64,
}

/// Bond parameters of one bar (all strictly validated by the caller; see [`BondSpec::validate`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondSpec {
    /// Total bar perimeter `n pi D` (mm).
    pub perimeter: f64,
    /// Design bond strength (MPa).
    pub fbd: f64,
    /// Elastic bond stiffness `kG Ec / D` (MPa per mm of slip).
    pub gb: f64,
    /// Plastic slope of the bond law as a fraction of `gb` (book default `1e-5`).
    pub rf: f64,
    /// Stop slip `alpha delta_u_max` (mm) for the bond and anchorage springs.
    pub slip_limit: f64,
    /// Turn angle (radians) above which a vertex is slip-locked.
    pub lock_angle: f64,
    /// Anchorage springs at the start and the end of the bar (`None` = a straight end, beta 0).
    pub ends: [Option<EndSpring>; 2],
}

impl BondSpec {
    /// `Err(reason)` unless every number is finite and in range.
    pub fn validate(&self) -> Result<(), &'static str> {
        let positive = |v: f64| v.is_finite() && v > 0.0;
        if !(positive(self.perimeter)
            && positive(self.fbd)
            && positive(self.gb)
            && positive(self.slip_limit))
        {
            return Err("bond perimeter, fbd, gb and slipLimit must be finite and positive");
        }
        if !(self.rf.is_finite() && (0.0..=1.0).contains(&self.rf)) {
            return Err("bond rf must lie in [0, 1]");
        }
        if !(self.lock_angle.is_finite() && self.lock_angle > 0.0) {
            return Err("bond lockAngle must be finite and positive");
        }
        for end in self.ends.iter().flatten() {
            if !(positive(end.fau)
                && positive(end.ka)
                && end.ksn.is_finite()
                && end.ksn > 0.0
                && end.ksn <= 1.0)
            {
                return Err("an anchorage end needs positive fau and ka and ksn in (0, 1]");
            }
        }
        Ok(())
    }
}

/// What a bonded bar adds to the model: the rod segments (to be added as `EmbeddedBar` elements, in order), the bond and
/// end springs, and where each bar point went.
#[derive(Debug)]
pub struct BondedBar {
    pub rods: Vec<EmbeddedBar>,
    /// Bond and end springs, to be added as `ZeroLength` elements.
    pub springs: Vec<BondSpring>,
    pub points: Vec<BondPoint>,
}

#[derive(Debug)]
pub struct BondSpring {
    pub element: ZeroLength,
    /// Index into [`BondedBar::points`].
    pub point: usize,
    pub end: bool,
    /// Tributary length (mm) for a bond spring; `0` for an end spring.
    pub tributary: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct BondPoint {
    pub coords: [f64; 2],
    pub image: NodeId,
    /// `None` at a slip-locked vertex (the segments attach to the image node).
    pub rod: Option<NodeId>,
    pub tangent: [f64; 2],
}

fn host_nodes(host: &EmbeddedBarHost) -> Vec<NodeId> {
    match host {
        EmbeddedBarHost::Tri3 { nodes, .. } => nodes.to_vec(),
        EmbeddedBarHost::Quad4 { nodes, .. } => nodes.to_vec(),
        EmbeddedBarHost::Rod { nodes, .. } => nodes.to_vec(),
    }
}

/// Builds the bonded model of one bar from its host-bound `segments` (as returned by `embed_polyline_in`, in path order).
/// Adds the image and rod nodes and the constraints to `domain`; returns the elements for the caller to register.
///
/// # Errors
///
/// [`ModelError::EmbeddedBarDegeneratePath`] with no segments. The caller validates `spec` (see [`BondSpec::validate`]).
pub fn bond_polyline(
    domain: &mut Domain,
    segments: &[EmbeddedBar],
    area: f64,
    material: &Material,
    spec: &BondSpec,
) -> Result<BondedBar, ModelError> {
    if segments.is_empty() {
        return Err(ModelError::EmbeddedBarDegeneratePath);
    }
    let count = segments.len() + 1;
    let mut coords = Vec::with_capacity(count);
    for segment in segments {
        let [start, _] = segment
            .end_points()
            .expect("embed_polyline sets end points");
        coords.push(start);
    }
    coords.push(
        segments[segments.len() - 1]
            .end_points()
            .expect("end points")[1],
    );

    // Tangent and turn angle at each point.
    let directions: Vec<[f64; 2]> = segments.iter().map(|s| s.direction()).collect();
    let mut tangents = Vec::with_capacity(count);
    let mut locked = vec![false; count];
    for k in 0..count {
        if k == 0 {
            tangents.push(directions[0]);
        } else if k == count - 1 {
            tangents.push(directions[count - 2]);
        } else {
            let (a, b) = (directions[k - 1], directions[k]);
            let sum = [a[0] + b[0], a[1] + b[1]];
            let norm = sum[0].hypot(sum[1]);
            // Opposite directions (a hairpin) have no tangent: lock the vertex, the bar doubles back on itself.
            let angle = (a[0] * b[0] + a[1] * b[1]).clamp(-1.0, 1.0).acos();
            locked[k] = angle > spec.lock_angle || norm < 1e-9;
            tangents.push(if norm < 1e-9 {
                a
            } else {
                [sum[0] / norm, sum[1] / norm]
            });
        }
    }

    // Image nodes tied to their host's shape functions.
    let mut points = Vec::with_capacity(count);
    for k in 0..count {
        let image = domain.add_node(Node::new(coords[k]).fix(2));
        let (segment, at_end) = if k < segments.len() {
            (&segments[k], false)
        } else {
            (&segments[k - 1], true)
        };
        let shapes = segment
            .end_shapes()
            .ok_or(ModelError::EmbeddedBarDegeneratePath)?;
        let values = if at_end { &shapes.1 } else { &shapes.0 };
        let hosts = host_nodes(segment.host());
        for dof in 0..2 {
            let terms: Vec<(NodeId, usize, f64)> = hosts
                .iter()
                .zip(values.iter())
                .map(|(&n, &v)| (n, dof, v))
                .collect();
            domain.add_constraint((image, dof), &terms);
        }
        let rod = if locked[k] {
            None
        } else {
            let rod = domain.add_node(Node::new(coords[k]));
            // Normal tie, solved for the DOF with the larger normal component: n . (u_R - u_C) = 0.
            let t = tangents[k];
            let n = [-t[1], t[0]];
            let (slave, other, ratio) = if n[1].abs() >= n[0].abs() {
                (1, 0, n[0] / n[1])
            } else {
                (0, 1, n[1] / n[0])
            };
            domain.add_constraint(
                (rod, slave),
                &[
                    (image, slave, 1.0),
                    (rod, other, -ratio),
                    (image, other, ratio),
                ],
            );
            Some(rod)
        };
        points.push(BondPoint {
            coords: coords[k],
            image,
            rod,
            tangent: tangents[k],
        });
    }

    let node_of = |k: usize| points[k].rod.unwrap_or(points[k].image);
    let mut rods = Vec::with_capacity(segments.len());
    for (s, segment) in segments.iter().enumerate() {
        let bar = EmbeddedBar::between_nodes(
            [node_of(s), node_of(s + 1)],
            directions[s],
            segment.length(),
            area,
            material.clone(),
        )
        .with_end_points(coords[s], coords[s + 1]);
        rods.push(bar);
    }

    let mut springs = Vec::new();
    for (k, point) in points.iter().enumerate() {
        let Some(rod) = point.rod else { continue };
        let neighbours = [k.checked_sub(1), (k < segments.len()).then_some(k)];
        let tributary: f64 = neighbours
            .iter()
            .flatten()
            .map(|&s| segments[s].length() / 2.0)
            .sum();
        let force_scale = spec.perimeter * tributary;
        let orientation =
            Orientation::in_plane(point.tangent[0], point.tangent[1]).expect("a unit tangent");
        let bond = ZeroLength::new(point.image, rod)
            .with_material(
                0,
                Material::steel01(
                    spec.fbd * force_scale,
                    spec.gb * force_scale,
                    spec.rf,
                    0.0,
                    1.0,
                    0.0,
                    1.0,
                ),
            )
            .with_slip_limit(0, spec.slip_limit, LimitKind::BondSlip)
            .with_orientation(orientation)
            .expect("an in-plane orientation");
        springs.push(BondSpring {
            element: bond,
            point: k,
            end: false,
            tributary,
        });
        let end = if k == 0 {
            spec.ends[0]
        } else if k == count - 1 {
            spec.ends[1]
        } else {
            None
        };
        if let Some(end) = end {
            let anchorage = ZeroLength::new(point.image, rod)
                .with_material(
                    0,
                    Material::steel01(end.fau, end.ka, end.ksn, 0.0, 1.0, 0.0, 1.0),
                )
                .with_slip_limit(0, spec.slip_limit, LimitKind::AnchorageSlip)
                .with_orientation(orientation)
                .expect("an in-plane orientation");
            springs.push(BondSpring {
                element: anchorage,
                point: k,
                end: true,
                tributary: 0.0,
            });
        }
    }
    Ok(BondedBar {
        rods,
        springs,
        points,
    })
}
