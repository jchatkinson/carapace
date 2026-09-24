//! Wire-format material arena and its resolution into `core::Material`
//! values. See pysees-handoff.md's "Materials/sections arena": composite
//! materials reference other arena entries by index rather than duplicating
//! them inline, mirroring `core`'s own arena-of-indices philosophy.

use carapace_core::model::Material;

use super::error::DecodeError;

/// One arena entry. Composite variants reference other entries by their
/// `u32` index into the same arena, resolved by [`resolve_materials`].
///
/// Deliberately narrower than `core::Material`'s full catalog for M10's
/// first decoder: `Hysteretic`/`Pinching4` are not yet supported (their
/// wire representation would need to mirror their large field lists, and
/// neither is needed by the first acceptance case).
#[derive(Debug, Clone)]
pub enum MaterialSpec {
    Elastic {
        e: f64,
    },
    ElasticPp {
        e: f64,
        eyp: f64,
    },
    Gap {
        e: f64,
        gap: f64,
    },
    Ent {
        e: f64,
    },
    Steel01 {
        fy: f64,
        e0: f64,
        b: f64,
        a1: f64,
        a2: f64,
        a3: f64,
        a4: f64,
    },
    Concrete01 {
        fpc: f64,
        epsc0: f64,
        fpcu: f64,
        epscu: f64,
    },
    #[allow(clippy::too_many_arguments)]
    Steel02 {
        fy: f64,
        e0: f64,
        b: f64,
        r0: f64,
        cr1: f64,
        cr2: f64,
        a1: f64,
        a2: f64,
        a3: f64,
        a4: f64,
    },
    Concrete02 {
        fc: f64,
        epsc0: f64,
        fcu: f64,
        epscu: f64,
        rat: f64,
        ft: f64,
        ets: f64,
    },
    /// Each child arena index paired with its parallel factor.
    Parallel {
        children: Vec<(u32, f64)>,
    },
    /// Child arena indices, combined with `core`'s default series
    /// tolerance/iteration cap (`Material::series`'s `tol=1e-10,
    /// max_iter=30`) — not yet wire-configurable.
    Series {
        children: Vec<u32>,
    },
    MinMax {
        inner: u32,
        min_strain: f64,
        max_strain: f64,
    },
}

/// Resolves every entry of a wire-format material arena into a `core`
/// `Material`, in arena order. Composite entries may reference any other
/// arena index (not just earlier ones); a reference cycle is rejected
/// explicitly rather than overflowing the call stack.
pub fn resolve_materials(arena: &[MaterialSpec]) -> Result<Vec<Material>, DecodeError> {
    let mut resolved: Vec<Option<Material>> = vec![None; arena.len()];
    for index in 0..arena.len() {
        resolve_one(arena, index, &mut resolved, &mut Vec::new())?;
    }
    Ok(resolved
        .into_iter()
        .map(|m| m.expect("every arena index visited above"))
        .collect())
}

fn resolve_one(
    arena: &[MaterialSpec],
    index: usize,
    resolved: &mut [Option<Material>],
    stack: &mut Vec<usize>,
) -> Result<Material, DecodeError> {
    if let Some(material) = &resolved[index] {
        return Ok(material.clone());
    }
    if stack.contains(&index) {
        return Err(DecodeError::CyclicMaterialReference {
            index: index as u32,
        });
    }
    let spec = arena.get(index).ok_or(DecodeError::UnknownMaterialIndex {
        table: "materials",
        row: index as u32,
    })?;

    stack.push(index);
    let material = match spec {
        MaterialSpec::Elastic { e } => Material::Elastic { e: *e },
        MaterialSpec::ElasticPp { e, eyp } => Material::elastic_pp(*e, *eyp),
        MaterialSpec::Gap { e, gap } => Material::Gap { e: *e, gap: *gap },
        MaterialSpec::Ent { e } => Material::Ent { e: *e },
        MaterialSpec::Steel01 {
            fy,
            e0,
            b,
            a1,
            a2,
            a3,
            a4,
        } => Material::steel01(*fy, *e0, *b, *a1, *a2, *a3, *a4),
        MaterialSpec::Concrete01 {
            fpc,
            epsc0,
            fpcu,
            epscu,
        } => Material::concrete01(*fpc, *epsc0, *fpcu, *epscu),
        MaterialSpec::Steel02 {
            fy,
            e0,
            b,
            r0,
            cr1,
            cr2,
            a1,
            a2,
            a3,
            a4,
        } => Material::steel02(*fy, *e0, *b, *r0, *cr1, *cr2, *a1, *a2, *a3, *a4),
        MaterialSpec::Concrete02 {
            fc,
            epsc0,
            fcu,
            epscu,
            rat,
            ft,
            ets,
        } => Material::concrete02(*fc, *epsc0, *fcu, *epscu, *rat, *ft, *ets),
        MaterialSpec::Parallel { children } => {
            let mut resolved_children = Vec::with_capacity(children.len());
            for (child_index, factor) in children {
                let child = resolve_one(arena, *child_index as usize, resolved, stack)?;
                resolved_children.push((child, *factor));
            }
            Material::Parallel(resolved_children)
        }
        MaterialSpec::Series { children } => {
            let resolved_children = children
                .iter()
                .map(|child_index| resolve_one(arena, *child_index as usize, resolved, stack))
                .collect::<Result<Vec<_>, _>>()?;
            Material::series(resolved_children)
        }
        MaterialSpec::MinMax {
            inner,
            min_strain,
            max_strain,
        } => {
            let inner = resolve_one(arena, *inner as usize, resolved, stack)?;
            Material::min_max(inner, *min_strain, *max_strain)
        }
    };
    stack.pop();

    resolved[index] = Some(material.clone());
    Ok(material)
}
