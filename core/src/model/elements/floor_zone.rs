//! Bay math of a stud wall's floor zone (`StudWallLayout::floor`): the depth between the framing top and the
//! exterior top row, as constant-shear bays plus an axial rim bar. Rows are `[bottom i, bottom j, top i, top j]`,
//! each with `(u, v)`, in wall-local axes; the bottom row is the framing top (`y = H`), the top row the exterior
//! nodes (`y = H + d`).

use nalgebra::{SMatrix, SVector};

/// Shear strain operator `gamma = (u_Ti + u_Tj - u_Bi - u_Bj)/(2d) + (v_Bj - v_Bi)/a`. The rotation term uses the
/// bottom row only, so the vertical couple `q a d` of the bay lands on the framing top (where the rim bears), and
/// a link across the zone carries the moment at the exterior top: the base of the wall above, as in O86 Fig. 11.7
/// (`R21 = F21 (H + h2)/L1`).
pub(crate) fn bay_operator(a: f64, d: f64) -> SVector<f64, 8> {
    let local = [
        (-1. / (2. * d), -1. / a),
        (-1. / (2. * d), 1. / a),
        (1. / (2. * d), 0.),
        (1. / (2. * d), 0.),
    ];
    SVector::from_fn(|j, _| {
        let (u, v) = local[j / 2];
        if j % 2 == 0 {
            u
        } else {
            v
        }
    })
}

/// Bay stiffness: shear energy `ga a d gamma²/2` (`ga = 0` for a bay with no shear field) plus the rim bar
/// between the two top nodes, `EA/a`.
pub(crate) fn bay_stiffness(a: f64, d: f64, ga: f64, rim_ea: f64) -> SMatrix<f64, 8, 8> {
    let b = bay_operator(a, d);
    let mut k = b * b.transpose() * (ga * a * d);
    let c = rim_ea / a;
    let q = [-1., 0., 1., 0.];
    for m in 0..4 {
        for n in 0..4 {
            k[(4 + m, 4 + n)] += c * q[m] * q[n];
        }
    }
    k
}

/// Shear strain, shear flow and rim force (tension positive) of a bay for local displacements `u`.
pub(crate) fn bay_response(a: f64, d: f64, ga: f64, rim_ea: f64, u: &[f64; 8]) -> (f64, f64, f64) {
    let gamma = bay_operator(a, d).dot(&SVector::from(*u));
    (gamma, ga * gamma, rim_ea / a * (u[6] - u[4]))
}
