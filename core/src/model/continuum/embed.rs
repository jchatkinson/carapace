//! `embed_polyline`: discretizes a rebar polyline into `EmbeddedBar`
//! segments against a 2D `Domain`'s continuum (`Tri3`/`Quad4`) elements
//! (plan `docs/csfm-plan.md` §2.3). A core model-building helper, not JS:
//! the wire format hands Carapace polylines, not precomputed shape
//! functions, and this is where the discretization itself is owned and
//! tested.
//!
//! Each drawn span is split at every mesh-edge crossing, so every
//! resulting segment lies in exactly one host; both of a segment's ends
//! are bound to the host containing the segment's own **midpoint**
//! (unlike CalcsApp's `rod.ts`, which binds each end point on its own —
//! identical on a conforming mesh, since neighbors share the edge's
//! nodes and off-edge shape functions vanish there, but this keeps every
//! segment within one host's `<= 4` nodes without a separate rule for an
//! end that sits exactly on a shared edge).

use super::{quad4_contains, quad4_natural_coords, tri3_contains, tri3_natural_coords};
use crate::model::{Domain, Element, EmbeddedBar, Material, ModelError};

const EPSILON: f64 = 1e-9;
/// Fixed tolerance for the natural-coordinate "is this point inside the
/// reference cell" test — natural coordinates are already dimensionless
/// (`[-1, 1]` or barycentric), so, unlike the physical tolerance below, no
/// per-host scaling is needed.
const NATURAL_TOL: f64 = 1e-7;
/// How far inside the reference cell a point must be to count as
/// "strictly interior" for ambiguous-host detection (see
/// `Host::contains_strict`) — comfortably above floating-point noise on a
/// boundary point, comfortably below any real element's extent.
const STRICT_MARGIN: f64 = 1e-6;
const MAX_NEWTON_ITERATIONS: usize = 50;

fn cross(ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    ax * by - ay * bx
}

/// The parameter `t` along `a -> b` where segment `a-b` properly crosses
/// segment `c-d` (both `t` and the crossing's own parameter along `c-d`
/// within `[0, 1]`, with `EPSILON` slack), or `None` if they are parallel
/// or do not cross. Ported from CalcsApp's `rod.ts::segmentParameters`.
fn segment_crossing(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> Option<f64> {
    let (rx, ry) = (b[0] - a[0], b[1] - a[1]);
    let (sx, sy) = (d[0] - c[0], d[1] - c[1]);
    let denominator = cross(rx, ry, sx, sy);
    if denominator.abs() < EPSILON {
        return None;
    }
    let (qx, qy) = (c[0] - a[0], c[1] - a[1]);
    let t = cross(qx, qy, sx, sy) / denominator;
    let u = cross(qx, qy, rx, ry) / denominator;
    if (-EPSILON..=1.0 + EPSILON).contains(&t) && (-EPSILON..=1.0 + EPSILON).contains(&u) {
        Some(t.clamp(0.0, 1.0))
    } else {
        None
    }
}

/// Every continuum element's edges, as node-coordinate pairs, in no
/// particular order and with duplicates for a shared edge (harmless here:
/// `unique_parameters` below already deduplicates the crossing parameters
/// those duplicates would produce).
fn host_edges(domain: &Domain) -> Vec<([f64; 2], [f64; 2])> {
    let mut edges = Vec::new();
    for (_, element) in domain.elements() {
        let coords: Vec<[f64; 2]> = match element {
            Element::Tri3(t) => t.nodes.iter().map(|&id| domain.node(id).coords).collect(),
            Element::Quad4(q) => q.nodes.iter().map(|&id| domain.node(id).coords).collect(),
            _ => continue,
        };
        for i in 0..coords.len() {
            edges.push((coords[i], coords[(i + 1) % coords.len()]));
        }
    }
    edges
}

/// Sorts and removes near-duplicate crossing parameters.
fn unique_parameters(mut values: Vec<f64>) -> Vec<f64> {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut unique: Vec<f64> = Vec::with_capacity(values.len());
    for value in values.drain(..) {
        if unique
            .last()
            .is_none_or(|&last| (value - last).abs() > EPSILON)
        {
            unique.push(value);
        }
    }
    unique
}

/// Splits `polyline` into the endpoints of every sub-span that lies
/// entirely in one mesh region (plan §2.3, §4 D): consecutive drawn
/// vertices, plus every point where the span crosses a host edge. Drops a
/// consecutive pair that is already coincident. Errs if fewer than two
/// vertices are given, or if every span was zero length.
fn discretize_polyline(
    edges: &[([f64; 2], [f64; 2])],
    polyline: &[[f64; 2]],
) -> Result<Vec<[f64; 2]>, ModelError> {
    if polyline.len() < 2 {
        return Err(ModelError::EmbeddedBarDegeneratePath);
    }
    let mut points: Vec<[f64; 2]> = Vec::new();
    for span in 0..polyline.len() - 1 {
        let (a, b) = (polyline[span], polyline[span + 1]);
        if (b[0] - a[0]).hypot(b[1] - a[1]) < EPSILON {
            continue;
        }
        let mut parameters = vec![0.0, 1.0];
        for &(c, d) in edges {
            if let Some(t) = segment_crossing(a, b, c, d) {
                parameters.push(t);
            }
        }
        let crossings = unique_parameters(parameters);
        for window in 0..crossings.len() - 1 {
            let (first, last) = (crossings[window], crossings[window + 1]);
            // No special-case skip for "the first point of a later span
            // duplicates the previous span's last point": `is_new` below
            // already de-duplicates against whatever was *actually* last
            // pushed, which is correct even when an earlier span
            // contributed no points at all (e.g. a zero-length span was
            // skipped above).
            for t in [first, last] {
                let point = [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];
                let is_new = points
                    .last()
                    .map(|last| (last[0] - point[0]).hypot(last[1] - point[1]) > EPSILON)
                    .unwrap_or(true);
                if is_new {
                    points.push(point);
                }
            }
        }
    }
    if points.len() < 2 {
        return Err(ModelError::EmbeddedBarDegeneratePath);
    }
    Ok(points)
}

/// One continuum host's node coordinates, found by the segment's natural-
/// coordinate type (`Tri3`/`Quad4`) and node list.
enum Host {
    Tri3 {
        nodes: [crate::model::NodeId; 3],
        coords: [[f64; 2]; 3],
    },
    Quad4 {
        nodes: [crate::model::NodeId; 4],
        coords: [[f64; 2]; 4],
    },
}

impl Host {
    fn natural_coords(&self, point: [f64; 2], physical_tol: f64) -> Option<[f64; 2]> {
        match self {
            Host::Tri3 { coords, .. } => tri3_natural_coords(coords, point).ok(),
            Host::Quad4 { coords, .. } => {
                quad4_natural_coords(coords, point, physical_tol, MAX_NEWTON_ITERATIONS).ok()
            }
        }
    }

    fn contains(&self, natural: [f64; 2]) -> bool {
        match self {
            Host::Tri3 { .. } => tri3_contains(natural, NATURAL_TOL),
            Host::Quad4 { .. } => quad4_contains(natural, NATURAL_TOL),
        }
    }

    /// Strictly interior (a negative tolerance shrinks the reference cell
    /// inward by `STRICT_MARGIN`): true only when the point is clearly
    /// inside, not merely within `NATURAL_TOL` of the boundary. Used to
    /// tell a genuinely overlapping host apart from two hosts that simply
    /// share an edge or a node — on a conforming mesh, a point exactly on
    /// a shared edge is a boundary match for *both* neighbors (that is
    /// fine: their shape functions agree there), but a point strictly
    /// interior to more than one host is a real ambiguity.
    fn contains_strict(&self, natural: [f64; 2]) -> bool {
        match self {
            Host::Tri3 { .. } => tri3_contains(natural, -STRICT_MARGIN),
            Host::Quad4 { .. } => quad4_contains(natural, -STRICT_MARGIN),
        }
    }

    fn shape(&self, natural: [f64; 2]) -> Vec<f64> {
        match self {
            Host::Tri3 { .. } => super::tri3_shape(natural[0], natural[1]).to_vec(),
            Host::Quad4 { .. } => super::quad4_shape(natural[0], natural[1]).to_vec(),
        }
    }

    /// A length scale for this host's own physical tolerance (plan §2.3:
    /// "separate tolerances: a physical one scaled to element size").
    fn scale(&self) -> f64 {
        let coords: &[[f64; 2]] = match self {
            Host::Tri3 { coords, .. } => coords,
            Host::Quad4 { coords, .. } => coords,
        };
        let (mut max_x, mut min_x, mut max_y, mut min_y) = (f64::MIN, f64::MAX, f64::MIN, f64::MAX);
        for &[x, y] in coords {
            max_x = max_x.max(x);
            min_x = min_x.min(x);
            max_y = max_y.max(y);
            min_y = min_y.min(y);
        }
        (max_x - min_x).max(max_y - min_y).max(1.0)
    }

    fn bounding_box_contains(&self, point: [f64; 2], tolerance: f64) -> bool {
        let coords: &[[f64; 2]] = match self {
            Host::Tri3 { coords, .. } => coords,
            Host::Quad4 { coords, .. } => coords,
        };
        let (mut max_x, mut min_x, mut max_y, mut min_y) = (f64::MIN, f64::MAX, f64::MIN, f64::MAX);
        for &[x, y] in coords {
            max_x = max_x.max(x);
            min_x = min_x.min(x);
            max_y = max_y.max(y);
            min_y = min_y.min(y);
        }
        point[0] >= min_x - tolerance
            && point[0] <= max_x + tolerance
            && point[1] >= min_y - tolerance
            && point[1] <= max_y + tolerance
    }
}

fn hosts(domain: &Domain) -> Vec<Host> {
    let mut found = Vec::new();
    for (_, element) in domain.elements() {
        match element {
            Element::Tri3(t) => found.push(Host::Tri3 {
                nodes: t.nodes,
                coords: std::array::from_fn(|i| domain.node(t.nodes[i]).coords),
            }),
            Element::Quad4(q) => found.push(Host::Quad4 {
                nodes: q.nodes,
                coords: std::array::from_fn(|i| domain.node(q.nodes[i]).coords),
            }),
            _ => {}
        }
    }
    found
}

/// Every host (brute force, bounding-box filtered — plan §2.3 says to
/// start there and add a spatial index only if benchmarks call for it)
/// whose natural-coordinate inverse of `point` lands inside its own
/// reference cell.
fn hosts_containing(hosts: &[Host], point: [f64; 2]) -> Vec<(&Host, [f64; 2], bool)> {
    let mut matches = Vec::new();
    for host in hosts {
        let tolerance = EPSILON * host.scale();
        if !host.bounding_box_contains(point, tolerance) {
            continue;
        }
        if let Some(natural) = host.natural_coords(point, tolerance) {
            if host.contains(natural) {
                matches.push((host, natural, host.contains_strict(natural)));
            }
        }
    }
    matches
}

/// Discretizes `polyline` (mm) into `EmbeddedBar`s of `area` (mm^2) and
/// `material`, against `domain`'s `Tri3`/`Quad4` elements.
///
/// # Errors
///
/// [`ModelError::EmbeddedBarDegeneratePath`] if `polyline` has fewer than
/// two vertices or no non-zero-length span; [`ModelError::
/// EmbeddedBarOutsideMesh`] if a resulting segment's midpoint lies inside
/// no host, or its own or its endpoints' natural-coordinate inverse does
/// not converge; [`ModelError::EmbeddedBarAmbiguousHost`] if a segment's
/// midpoint lies inside more than one host.
pub fn embed_polyline(
    domain: &Domain,
    polyline: &[[f64; 2]],
    area: f64,
    material: &Material,
) -> Result<Vec<EmbeddedBar>, ModelError> {
    let edges = host_edges(domain);
    let points = discretize_polyline(&edges, polyline)?;
    let hosts = hosts(domain);

    let mut bars = Vec::with_capacity(points.len() - 1);
    for segment in 0..points.len() - 1 {
        let (a, b) = (points[segment], points[segment + 1]);
        let midpoint = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];

        let matches = hosts_containing(&hosts, midpoint);
        if matches.is_empty() {
            return Err(ModelError::EmbeddedBarOutsideMesh { segment });
        }
        let strictly_interior: Vec<_> = matches.iter().filter(|&&(_, _, strict)| strict).collect();
        if strictly_interior.len() > 1 {
            return Err(ModelError::EmbeddedBarAmbiguousHost { segment });
        }
        // Prefer a strictly interior match; otherwise (every match is a
        // boundary coincidence between neighbors) any one of them gives
        // the same physically meaningful result.
        let (host, _, _) = *strictly_interior.first().copied().unwrap_or(&matches[0]);

        let tolerance = EPSILON * host.scale();
        let natural_a = host
            .natural_coords(a, tolerance)
            .filter(|&n| host.contains(n))
            .ok_or(ModelError::EmbeddedBarOutsideMesh { segment })?;
        let natural_b = host
            .natural_coords(b, tolerance)
            .filter(|&n| host.contains(n))
            .ok_or(ModelError::EmbeddedBarOutsideMesh { segment })?;

        let n_a = host.shape(natural_a);
        let n_b = host.shape(natural_b);
        let length = (b[0] - a[0]).hypot(b[1] - a[1]);
        let direction = [(b[0] - a[0]) / length, (b[1] - a[1]) / length];

        bars.push(match host {
            Host::Tri3 { nodes, .. } => EmbeddedBar::in_tri3(
                *nodes,
                std::array::from_fn(|i| n_a[i]),
                std::array::from_fn(|i| n_b[i]),
                direction,
                length,
                area,
                material.clone(),
            ),
            Host::Quad4 { nodes, .. } => EmbeddedBar::in_quad4(
                *nodes,
                std::array::from_fn(|i| n_a[i]),
                std::array::from_fn(|i| n_b[i]),
                direction,
                length,
                area,
                material.clone(),
            ),
        });
    }
    Ok(bars)
}
