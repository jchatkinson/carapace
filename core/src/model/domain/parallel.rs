//! Threaded element loop for the pattern assembly (`parallel` feature, native only).
//!
//! Workers evaluate elements in chunks, each into its own flat arena of recorded contributions (no
//! per-element allocation); the arenas are then replayed into the real sink on the calling thread in element
//! order, with the same arithmetic as the serial `PatternSink::add`. The scatter order is unchanged, so the
//! assembled matrix and resistance are bit-identical to the serial loop for any thread count.

use crate::model::NodeKey;
use nalgebra::{SMatrix, SVector};
use rayon::prelude::*;
use slotmap::KeyData;

use super::*;

/// Below this many elements the fork/join overhead outweighs the saving.
const MIN_ELEMENTS: usize = 256;
const CHUNK: usize = 128;

/// One chunk's recorded contributions: `n` dofs (as `(node ffi key, slot)`), then `n*n` stiffness entries
/// row-major, then `n` resistance entries, per contribution.
#[derive(Default)]
struct Arena {
    sizes: Vec<usize>,
    dofs: Vec<(u64, u8)>,
    values: Vec<f64>,
}

struct RecordSink<'a, const NDIM: usize, const NDOF: usize, NId, E>
where
    NId: NodeKey,
    E: ElementOps<NDIM, NDOF, NId>,
{
    #[cfg_attr(not(debug_assertions), allow(dead_code))]
    domain: &'a Domain<NDIM, NDOF, NId, E>,
    arena: Arena,
}

impl<const NDIM: usize, const NDOF: usize, NId, E> TangentSink<NId>
    for RecordSink<'_, NDIM, NDOF, NId, E>
where
    NId: NodeKey,
    E: ElementOps<NDIM, NDOF, NId>,
{
    fn add<const N: usize>(
        &mut self,
        dofs: &[DofRef<NId>; N],
        k: &SMatrix<f64, N, N>,
        r: &SVector<f64, N>,
    ) {
        #[cfg(debug_assertions)]
        self.domain.assert_declared(dofs, Some(k), r.as_slice());
        let arena = &mut self.arena;
        arena.sizes.push(N);
        arena.dofs.extend(
            dofs.iter()
                .map(|&(node, slot)| (node.data().as_ffi(), slot)),
        );
        for a in 0..N {
            for b in 0..N {
                arena.values.push(k[(a, b)]);
            }
        }
        arena.values.extend(r.iter().copied());
    }
}

/// Fills `sink` from the elements on the configured pool. `false` (nothing done) when the build is serial,
/// the model is small, or no pool is set.
pub(super) fn fill<const NDIM: usize, const NDOF: usize, NId, E>(
    domain: &Domain<NDIM, NDOF, NId, E>,
    active: &[(f64, &LoadPattern<NDOF, NId, E::Id, E::Load>)],
    view: &NodeView<'_, NDIM, NDOF, NId>,
    sink: &mut PatternSink<'_, NDIM, NDOF, NId, E>,
) -> bool
where
    NId: NodeKey,
    E: ElementOps<NDIM, NDOF, NId>,
{
    let Some(pool) = crate::analysis::pool() else {
        return false;
    };
    if domain.elements.len() < MIN_ELEMENTS {
        return false;
    }
    let elements: Vec<(E::Id, &E)> = domain.elements.iter().collect();
    let arenas: Vec<Arena> = pool.install(|| {
        elements
            .par_chunks(CHUNK)
            .map(|chunk| {
                let mut record = RecordSink {
                    domain,
                    arena: Arena::default(),
                };
                for &(id, element) in chunk {
                    let load = effective_element_load(active, id);
                    element.assemble_tangent(view, load.as_ref(), &mut record);
                }
                record.arena
            })
            .collect()
    });
    for arena in &arenas {
        let (mut at_dof, mut at_value) = (0, 0);
        for &n in &arena.sizes {
            let dofs: Vec<DofRef<NId>> = arena.dofs[at_dof..at_dof + n]
                .iter()
                .map(|&(key, slot)| (NId::from(KeyData::from_ffi(key)), slot))
                .collect();
            let k = &arena.values[at_value..at_value + n * n];
            let r = &arena.values[at_value + n * n..at_value + n * n + n];
            sink.add_recorded(&dofs, k, r);
            at_dof += n;
            at_value += n * n + n;
        }
    }
    true
}

impl<const NDIM: usize, const NDOF: usize, NId, E> PatternSink<'_, NDIM, NDOF, NId, E>
where
    NId: NodeKey,
    E: ElementOps<NDIM, NDOF, NId>,
{
    /// `TangentSink::add` for a recorded contribution (`k` row-major `n x n`): same loops and the same
    /// floating-point expressions, so the result matches the serial path exactly.
    fn add_recorded(&mut self, dofs: &[DofRef<NId>], k: &[f64], r: &[f64]) {
        let n = dofs.len();
        for (a, &(node_a, slot_a)) in dofs.iter().enumerate() {
            for &(eq_a, coeff_a) in self.domain.dof_terms(node_a, slot_a as usize).iter() {
                self.resistance[eq_a] += coeff_a * r[a];
                for (b, &(node_b, slot_b)) in dofs.iter().enumerate() {
                    for &(eq_b, coeff_b) in self.domain.dof_terms(node_b, slot_b as usize).iter() {
                        match self.pattern.position(eq_a, eq_b) {
                            Some(at) => self.values[at] += coeff_a * coeff_b * k[a * n + b],
                            None => self.complete = false,
                        }
                    }
                }
            }
        }
    }
}
