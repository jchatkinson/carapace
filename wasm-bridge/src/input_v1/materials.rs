//! Wire-format material arena and its resolution into `core::Material`
//! values. See pysees-handoff.md's "Materials/sections arena": composite
//! materials reference other arena entries by index rather than duplicating
//! them inline, mirroring `core`'s own arena-of-indices philosophy.

use carapace_core::model::{Material, Pinching4DmgCyc};
use serde::{Deserialize, Serialize};

use super::error::DecodeError;

/// Wire-format mirror of `core::Pinching4DmgCyc` (which carries no `serde`
/// derives of its own — `core` stays free of any wasm/serde awareness, see
/// `decode.rs`'s module doc comment).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum Pinching4DmgCycSpec {
    EnergyBased,
    CycleBased,
}

/// One arena entry. Composite variants reference other entries by their
/// `u32` index into the same arena, resolved by [`resolve_materials`].
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
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
    /// Trilinear moment-rotation envelope with pinching/stiffness/strength
    /// degradation — see `core::Material::hysteretic`'s doc comment. Field
    /// order matches the constructor's positional parameters.
    #[allow(clippy::too_many_arguments)]
    Hysteretic {
        mom1p: f64,
        rot1p: f64,
        mom2p: f64,
        rot2p: f64,
        mom3p: f64,
        rot3p: f64,
        mom1n: f64,
        rot1n: f64,
        mom2n: f64,
        rot2n: f64,
        mom3n: f64,
        rot3n: f64,
        pinch_x: f64,
        pinch_y: f64,
        damfc1: f64,
        damfc2: f64,
        beta: f64,
    },
    /// Quadrilinear backbone with pinching and cyclic damage — see
    /// `core::Material::pinching4`'s doc comment. Field order matches the
    /// constructor's positional parameters.
    #[allow(clippy::too_many_arguments)]
    Pinching4 {
        stress1p: f64,
        strain1p: f64,
        stress2p: f64,
        strain2p: f64,
        stress3p: f64,
        strain3p: f64,
        stress4p: f64,
        strain4p: f64,
        stress1n: f64,
        strain1n: f64,
        stress2n: f64,
        strain2n: f64,
        stress3n: f64,
        strain3n: f64,
        stress4n: f64,
        strain4n: f64,
        r_disp_p: f64,
        r_force_p: f64,
        u_force_p: f64,
        r_disp_n: f64,
        r_force_n: f64,
        u_force_n: f64,
        gamma_k_params: [f64; 4],
        gamma_k_limit: f64,
        gamma_d_params: [f64; 4],
        gamma_d_limit: f64,
        gamma_f_params: [f64; 4],
        gamma_f_limit: f64,
        gamma_e: f64,
        dmg_cyc: Pinching4DmgCycSpec,
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
        MaterialSpec::Hysteretic {
            mom1p,
            rot1p,
            mom2p,
            rot2p,
            mom3p,
            rot3p,
            mom1n,
            rot1n,
            mom2n,
            rot2n,
            mom3n,
            rot3n,
            pinch_x,
            pinch_y,
            damfc1,
            damfc2,
            beta,
        } => Material::hysteretic(
            *mom1p, *rot1p, *mom2p, *rot2p, *mom3p, *rot3p, *mom1n, *rot1n, *mom2n, *rot2n,
            *mom3n, *rot3n, *pinch_x, *pinch_y, *damfc1, *damfc2, *beta,
        ),
        MaterialSpec::Pinching4 {
            stress1p,
            strain1p,
            stress2p,
            strain2p,
            stress3p,
            strain3p,
            stress4p,
            strain4p,
            stress1n,
            strain1n,
            stress2n,
            strain2n,
            stress3n,
            strain3n,
            stress4n,
            strain4n,
            r_disp_p,
            r_force_p,
            u_force_p,
            r_disp_n,
            r_force_n,
            u_force_n,
            gamma_k_params,
            gamma_k_limit,
            gamma_d_params,
            gamma_d_limit,
            gamma_f_params,
            gamma_f_limit,
            gamma_e,
            dmg_cyc,
        } => Material::pinching4(
            *stress1p, *strain1p, *stress2p, *strain2p, *stress3p, *strain3p, *stress4p,
            *strain4p, *stress1n, *strain1n, *stress2n, *strain2n, *stress3n, *strain3n,
            *stress4n, *strain4n, *r_disp_p, *r_force_p, *u_force_p, *r_disp_n, *r_force_n,
            *u_force_n, *gamma_k_params, *gamma_k_limit, *gamma_d_params, *gamma_d_limit,
            *gamma_f_params, *gamma_f_limit, *gamma_e,
            match dmg_cyc {
                Pinching4DmgCycSpec::EnergyBased => Pinching4DmgCyc::EnergyBased,
                Pinching4DmgCycSpec::CycleBased => Pinching4DmgCyc::CycleBased,
            },
        ),
    };
    stack.pop();

    resolved[index] = Some(material.clone());
    Ok(material)
}
