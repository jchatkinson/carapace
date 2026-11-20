//! Shared helpers for driving a `DisplacementControl` cyclic protocol on a
//! single controlled DOF — used by `tests/analysis/cyclic_protocol.rs`
//! and `examples/cyclic_material.rs` (the openseespy comparison harness in
//! `comparison/`).
//!
//! Also the generic tangent-vs-central-difference checker
//! (`docs/csfm-plan.md` phase 0): for a `Material`, a `PlaneMaterial`, or
//! any `ElementOps` element, perturb the input, central-difference the
//! returned stress/resistance, and compare against the returned tangent.
//! This is a smooth-branch check only — every caller supplies a `kink`
//! predicate so the checker can refuse a stencil that straddles a regime
//! or yield kink, rather than silently averaging across one (the defect
//! `docs/csfm-plan.md` §1 and §4.A note TS's own finite-difference stencil
//! had). Numerical policies (seeded tangents at isolated points, clamped
//! tangents, rupture secants) are deliberate and are tested directly
//! against their specified value instead, never through this checker.

use nalgebra::{DMatrix, DVector, SMatrix, SVector};
use slotmap::{Key, SlotMap};

use crate::analysis::{Analysis, Integrator};
use crate::model::{
    DofRef, ElementOps, Material, Node, NodeId, NodeView, PlaneMaterial, PlaneMatrix, PlaneVector,
    TangentSink,
};

/// One `DisplacementControl` leg via `Analysis::set_integrator` (no
/// rebuild) targeting an *absolute* displacement — the increment passed to
/// `DisplacementControl` is computed fresh from wherever the controlled DOF
/// actually is right now, not assumed from the previous leg's target.
pub fn leg(analysis: &mut Analysis, node: NodeId, dof: usize, target: f64) -> f64 {
    let current = analysis.domain().node(node).displacement[dof];
    analysis.set_integrator(Integrator::DisplacementControl {
        node,
        dof,
        increment: target - current,
    });
    analysis
        .step()
        .expect("cyclic protocol leg should converge")
        .load_factor
}

/// Walks the controlled DOF from wherever it is to `target` in steps no
/// larger than `step`, via repeated `leg` calls. Calls `on_point` after
/// every leg with `(displacement, load_factor)` so callers can record a
/// full trace, not just the endpoint.
pub fn ramp_to(
    analysis: &mut Analysis,
    node: NodeId,
    dof: usize,
    target: f64,
    step: f64,
    mut on_point: impl FnMut(f64, f64),
) {
    loop {
        let current = analysis.domain().node(node).displacement[dof];
        let remaining = target - current;
        if remaining.abs() < 1e-12 {
            return;
        }
        let next = current + remaining.clamp(-step, step);
        let load_factor = leg(analysis, node, dof, next);
        on_point(next, load_factor);
    }
}

/// Runs a reverse-cyclic protocol: for each `(amplitude, cycles)` tier, in
/// order, ramps 0 -> +amplitude -> -amplitude, `cycles` times, then back to
/// 0 before the next tier. Calls `on_point` after every leg — this is the
/// full recorded trace a stress-strain (or force-displacement) plot is
/// built from.
pub fn run_cyclic_protocol(
    analysis: &mut Analysis,
    node: NodeId,
    dof: usize,
    tiers: &[(f64, usize)],
    step: f64,
    mut on_point: impl FnMut(f64, f64),
) {
    for &(amplitude, cycles) in tiers {
        for _ in 0..cycles {
            ramp_to(analysis, node, dof, amplitude, step, &mut on_point);
            ramp_to(analysis, node, dof, -amplitude, step, &mut on_point);
        }
        ramp_to(analysis, node, dof, 0.0, step, &mut on_point);
    }
}

// ------------------------------------------------------- tangent-vs-FD checker

/// Step size and tolerance for [`assert_material_tangent_matches_fd`] and
/// friends. The central-difference step scales with the magnitude of the
/// point being perturbed (`h = max(h_min, h_rel * x.abs())`), since a fixed
/// absolute step is either too coarse at large strains or swamped by
/// cancellation at small ones. Comparisons are componentwise absolute-plus-
/// relative, as `docs/csfm-plan.md` §4 specifies for every CSFM tolerance:
/// `|fd - analytic| <= tol_abs + tol_rel * max(|fd|, |analytic|)`.
#[derive(Debug, Clone, Copy)]
pub struct FdConfig {
    pub h_rel: f64,
    pub h_min: f64,
    pub tol_abs: f64,
    pub tol_rel: f64,
}

impl Default for FdConfig {
    /// `rel 1e-7`, matching the plan's default for a law's own smooth-branch
    /// tangent against this same checker (`docs/csfm-plan.md` §4.A).
    fn default() -> Self {
        FdConfig {
            h_rel: 1e-6,
            h_min: 1e-9,
            tol_abs: 1e-9,
            tol_rel: 1e-7,
        }
    }
}

impl FdConfig {
    fn step(&self, x: f64) -> f64 {
        self.h_min.max(self.h_rel * x.abs())
    }

    fn close(&self, fd: f64, analytic: f64) -> bool {
        (fd - analytic).abs() <= self.tol_abs + self.tol_rel * fd.abs().max(analytic.abs())
    }
}

/// Central-differences `f` at `x` with a step scaled to `x`'s magnitude
/// (`FdConfig::step`), refusing (returning `None`) if `kink` reports that
/// the center or either shifted point sits at/near a branch boundary — the
/// plan's "refuses a stencil that straddles a regime or yield kink".
fn central_difference(
    x: f64,
    config: &FdConfig,
    kink: impl Fn(f64) -> bool,
    f: impl Fn(f64) -> f64,
) -> Option<f64> {
    let h = config.step(x);
    if kink(x - h) || kink(x) || kink(x + h) {
        return None;
    }
    Some((f(x + h) - f(x - h)) / (2.0 * h))
}

/// Checks `material`'s tangent at `strain` against a central difference of
/// its own stress, provided the stencil (`strain` and its two shifted
/// points) stays clear of every point `kink` reports as a regime or yield
/// boundary. Returns `true` if the check ran and passed, `false` if `kink`
/// made the checker refuse the stencil (the caller then knows to pick
/// another point, or to test that boundary directly instead — see
/// `docs/csfm-plan.md` §4.A/B, "numerical policies, tested on their own").
///
/// `kink(x)` should answer for a *trial* strain `x` evaluated against the
/// same committed `material` — e.g. for `Steel01`/`Concrete01` in
/// `core/tests/materials`, a yield or envelope-corner strain.
pub fn assert_material_tangent_matches_fd(
    material: &Material,
    strain: f64,
    config: &FdConfig,
    kink: impl Fn(f64) -> bool,
) -> bool {
    let (_, analytic) = material.trial_stress_tangent(strain);
    let Some(fd) = central_difference(strain, config, kink, |x| material.trial_stress_tangent(x).0)
    else {
        return false;
    };
    assert!(
        config.close(fd, analytic),
        "tangent {analytic} does not match central-difference {fd} at strain {strain} (h = {})",
        config.step(strain)
    );
    true
}

/// The `PlaneMaterial` analogue of [`assert_material_tangent_matches_fd`]:
/// perturbs each of the three engineering strain components independently
/// and compares each resulting column of the central-difference Jacobian
/// against the matching column of `D`. `kink(strain)` answers for a trial
/// *strain vector*, since a plane law's branch structure (e.g. CSFM's
/// principal-strain regimes) is not a function of one component alone.
pub fn assert_plane_material_tangent_matches_fd(
    material: &PlaneMaterial,
    strain: &PlaneVector,
    config: &FdConfig,
    kink: impl Fn(&PlaneVector) -> bool,
) -> bool {
    if kink(strain) {
        return false;
    }
    let (_, analytic) = material.trial_stress_tangent(strain);
    let mut fd = PlaneMatrix::zeros();
    for a in 0..3 {
        let h = config.step(strain[a]);
        let mut plus = *strain;
        plus[a] += h;
        let mut minus = *strain;
        minus[a] -= h;
        if kink(&plus) || kink(&minus) {
            return false;
        }
        let d_plus = material.trial_stress_tangent(&plus).0;
        let d_minus = material.trial_stress_tangent(&minus).0;
        fd.set_column(a, &((d_plus - d_minus) / (2.0 * h)));
    }
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                config.close(fd[(i, j)], analytic[(i, j)]),
                "D[{i},{j}] = {} does not match central-difference {} at strain {strain:?}",
                analytic[(i, j)],
                fd[(i, j)]
            );
        }
    }
    true
}

/// Captures one `assemble_tangent` call's dofs, tangent and resistance as
/// dense (heap) matrices, regardless of the element's own fixed-size `N` —
/// the generic element checker below is written once for every element,
/// which have different `N`.
struct Collect<NId> {
    dofs: Vec<DofRef<NId>>,
    k: DMatrix<f64>,
    r: DVector<f64>,
}

impl<NId: Copy> Default for Collect<NId> {
    fn default() -> Self {
        Collect {
            dofs: Vec::new(),
            k: DMatrix::zeros(0, 0),
            r: DVector::zeros(0),
        }
    }
}

impl<NId: Copy> TangentSink<NId> for Collect<NId> {
    fn add<const N: usize>(
        &mut self,
        dofs: &[DofRef<NId>; N],
        k: &SMatrix<f64, N, N>,
        r: &SVector<f64, N>,
    ) {
        self.dofs = dofs.to_vec();
        self.k = DMatrix::from_fn(N, N, |i, j| k[(i, j)]);
        self.r = DVector::from_fn(N, |i, _| r[i]);
    }
}

/// Checks one `ElementOps` element's `assemble_tangent` at the nodal state
/// currently in `nodes`: perturbs each of the element's own DOFs in turn
/// and central-differences the resulting resisting-force vector against
/// the matching column of the tangent. Leaves `nodes` exactly as given
/// (every perturbation is undone before returning).
///
/// There is no `kink` hook here: an element's smoothness depends on its
/// material(s) along a path through strain space that the element itself
/// determines, which this generic helper cannot predict. Callers drive
/// elements at a nodal state they already know sits on a smooth branch of
/// every embedded material (as the plan's element tests do), and test a
/// branch boundary directly instead of through this checker.
pub fn assert_element_tangent_matches_fd<const NDIM: usize, const NDOF: usize, NId, E>(
    element: &mut E,
    nodes: &mut SlotMap<NId, Node<NDIM, NDOF>>,
    load: Option<&E::Load>,
    config: &FdConfig,
) where
    NId: Key,
    E: ElementOps<NDIM, NDOF, NId>,
{
    element.prepare(&NodeView::new(nodes));
    let assemble = |nodes: &SlotMap<NId, Node<NDIM, NDOF>>| -> Collect<NId> {
        let view = NodeView::new(nodes);
        let mut sink = Collect::default();
        element.assemble_tangent(&view, load, &mut sink);
        sink
    };
    let base = assemble(nodes);
    let n = base.dofs.len();
    assert_eq!(base.k.nrows(), n);
    for a in 0..n {
        let (node_id, slot) = base.dofs[a];
        let slot = slot as usize;
        let x0 = nodes[node_id].displacement[slot];
        let h = config.step(x0).max(config.h_min);
        nodes[node_id].displacement[slot] = x0 + h;
        let plus = assemble(nodes);
        nodes[node_id].displacement[slot] = x0 - h;
        let minus = assemble(nodes);
        nodes[node_id].displacement[slot] = x0;
        for i in 0..n {
            let fd = (plus.r[i] - minus.r[i]) / (2.0 * h);
            let analytic = base.k[(i, a)];
            assert!(
                config.close(fd, analytic),
                "K[{i},{a}] = {analytic} does not match central-difference {fd} of resistance (dof {a} = {node_id:?}.{slot}, h = {h})"
            );
        }
    }
}
