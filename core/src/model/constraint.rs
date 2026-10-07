//! General linear multi-point constraints: one slave DOF equals a linear
//! combination of master DOFs on any nodes, `u_slave = sum(coeff * u_master)`.
//! `equal_dof`, the rigid diaphragms and rigid links are all generated as
//! these. Constraints are homogeneous (no offset term), so an increment of a
//! slave follows its masters' increments exactly.

use std::collections::{HashMap, HashSet};

use slotmap::Key;

use super::dof_table::{DofEntry, DofTable};

#[derive(Debug, Clone)]
pub(crate) struct LinearConstraint<NId> {
    pub(crate) slave: (NId, usize),
    /// `(node, dof, coefficient)`, already normalized (see `normalize`).
    pub(crate) terms: Vec<(NId, usize, f64)>,
}

/// Merges repeated masters (summing their coefficients) and drops terms
/// that vanish, exactly or relative to the largest coefficient. The order
/// of first appearance is kept so resolved terms are reproducible.
pub(crate) fn normalize<NId: Key>(terms: &[(NId, usize, f64)]) -> Vec<(NId, usize, f64)> {
    let mut merged: Vec<(NId, usize, f64)> = Vec::with_capacity(terms.len());
    for &(node, dof, coeff) in terms {
        match merged.iter_mut().find(|(n, d, _)| *n == node && *d == dof) {
            Some(existing) => existing.2 += coeff,
            None => merged.push((node, dof, coeff)),
        }
    }
    let largest = merged.iter().map(|t| t.2.abs()).fold(0.0_f64, f64::max);
    merged.retain(|t| t.2 != 0.0 && t.2.abs() > 1e-14 * largest);
    merged
}

/// Why `resolve` could not finish.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ResolveError<NId> {
    /// Two constraints define the same slave DOF.
    DuplicateSlave((NId, usize)),
    /// A slave depends on itself, directly or through other slaves.
    Cycle((NId, usize)),
}

/// Resolves every slave DOF to terms over the free equations in `table`.
///
/// `table` must already hold `Free`/`Fixed`/`Inactive` for every DOF that is
/// not a slave. A master that is itself a slave is substituted (chains of any
/// depth); a fixed or inactive master contributes nothing. A slave that
/// resolves to exactly `1 * u_eq` becomes an alias of that equation
/// (`DofEntry::Free`), which is also what an identity tie always was; one that
/// resolves to nothing is `Fixed`; anything else is `Affine`.
pub(crate) fn resolve<NId: Key, const NDOF: usize>(
    constraints: &[LinearConstraint<NId>],
    table: &mut DofTable<NId, NDOF>,
) -> Result<(), ResolveError<NId>> {
    let mut index: HashMap<(NId, usize), usize> = HashMap::with_capacity(constraints.len());
    for (i, constraint) in constraints.iter().enumerate() {
        if index.insert(constraint.slave, i).is_some() {
            return Err(ResolveError::DuplicateSlave(constraint.slave));
        }
    }

    let mut memo: HashMap<(NId, usize), Vec<(usize, f64)>> = HashMap::new();
    let mut visiting: HashSet<(NId, usize)> = HashSet::new();
    for constraint in constraints {
        let terms = resolve_slave(
            constraint.slave,
            constraints,
            &index,
            table,
            &mut memo,
            &mut visiting,
        )?;
        let entry = match terms.as_slice() {
            [] => DofEntry::Fixed,
            [(eq, coeff)] if *coeff == 1.0 => DofEntry::Free(*eq),
            _ => table.push_terms(&terms),
        };
        table.set(constraint.slave.0, constraint.slave.1, entry);
    }
    Ok(())
}

fn resolve_slave<NId: Key, const NDOF: usize>(
    slave: (NId, usize),
    constraints: &[LinearConstraint<NId>],
    index: &HashMap<(NId, usize), usize>,
    table: &DofTable<NId, NDOF>,
    memo: &mut HashMap<(NId, usize), Vec<(usize, f64)>>,
    visiting: &mut HashSet<(NId, usize)>,
) -> Result<Vec<(usize, f64)>, ResolveError<NId>> {
    if let Some(done) = memo.get(&slave) {
        return Ok(done.clone());
    }
    if !visiting.insert(slave) {
        return Err(ResolveError::Cycle(slave));
    }
    let mut terms: Vec<(usize, f64)> = Vec::new();
    let add = |terms: &mut Vec<(usize, f64)>, eq: usize, coeff: f64| match terms
        .iter_mut()
        .find(|(e, _)| *e == eq)
    {
        Some(existing) => existing.1 += coeff,
        None => terms.push((eq, coeff)),
    };
    for &(node, dof, coeff) in &constraints[index[&slave]].terms {
        if index.contains_key(&(node, dof)) {
            let inner = resolve_slave((node, dof), constraints, index, table, memo, visiting)?;
            for (eq, inner_coeff) in inner {
                add(&mut terms, eq, coeff * inner_coeff);
            }
        } else if let DofEntry::Free(eq) = table.entry(node, dof) {
            add(&mut terms, eq, coeff);
        }
    }
    terms.retain(|t| t.1 != 0.0);
    visiting.remove(&slave);
    memo.insert(slave, terms.clone());
    Ok(terms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::NodeId;
    use slotmap::SlotMap;

    fn nodes(n: usize) -> (SlotMap<NodeId, ()>, Vec<NodeId>) {
        let mut map = SlotMap::with_key();
        let ids = (0..n).map(|_| map.insert(())).collect();
        (map, ids)
    }

    #[test]
    fn normalize_merges_repeated_masters_and_drops_vanishing_terms() {
        let (_, n) = nodes(2);
        let terms = normalize(&[
            (n[0], 0, 1.0),
            (n[1], 2, 0.5),
            (n[0], 0, 2.0),
            (n[1], 2, -0.5),
            (n[1], 1, 1e-20),
        ]);
        assert_eq!(terms, vec![(n[0], 0, 3.0)]);
    }

    #[test]
    fn chains_substitute_and_unit_single_terms_alias() {
        let (map, n) = nodes(4);
        let mut table = DofTable::<NodeId, 3>::for_nodes(map.keys());
        table.set(n[0], 0, DofEntry::Free(0));
        table.set(n[0], 2, DofEntry::Free(1));
        // n3 = 2*n2, n2 = n1, n1 = n0.ux - 3*n0.rz
        let constraints = vec![
            LinearConstraint {
                slave: (n[3], 0),
                terms: vec![(n[2], 0, 2.0)],
            },
            LinearConstraint {
                slave: (n[2], 0),
                terms: vec![(n[1], 0, 1.0)],
            },
            LinearConstraint {
                slave: (n[1], 0),
                terms: vec![(n[0], 0, 1.0), (n[0], 2, -3.0)],
            },
        ];
        resolve(&constraints, &mut table).expect("acyclic");
        let terms_of = |node: NodeId| match table.entry(node, 0) {
            DofEntry::Affine { start, len } => table.terms(start, len).to_vec(),
            other => panic!("expected affine, got {other:?}"),
        };
        assert_eq!(terms_of(n[1]), vec![(0, 1.0), (1, -3.0)]);
        assert_eq!(terms_of(n[2]), vec![(0, 1.0), (1, -3.0)]);
        assert_eq!(terms_of(n[3]), vec![(0, 2.0), (1, -6.0)]);

        // A unit single term is just an alias of the master equation.
        let mut table = DofTable::<NodeId, 3>::for_nodes(map.keys());
        table.set(n[0], 0, DofEntry::Free(0));
        let tie = vec![LinearConstraint {
            slave: (n[1], 0),
            terms: vec![(n[0], 0, 1.0)],
        }];
        resolve(&tie, &mut table).unwrap();
        assert_eq!(table.entry(n[1], 0), DofEntry::Free(0));
    }

    #[test]
    fn a_coefficient_other_than_one_is_not_an_alias_and_fixed_masters_vanish() {
        let (map, n) = nodes(3);
        let mut table = DofTable::<NodeId, 3>::for_nodes(map.keys());
        table.set(n[0], 0, DofEntry::Free(0));
        let constraints = vec![
            LinearConstraint {
                slave: (n[1], 0),
                terms: vec![(n[0], 0, 2.0)],
            },
            // n[2].uy is Fixed in the table, so this slave resolves to nothing.
            LinearConstraint {
                slave: (n[2], 0),
                terms: vec![(n[2], 1, 5.0)],
            },
        ];
        resolve(&constraints, &mut table).unwrap();
        let DofEntry::Affine { start, len } = table.entry(n[1], 0) else {
            panic!("a coefficient of 2 must stay affine");
        };
        assert_eq!(table.terms(start, len), &[(0, 2.0)]);
        assert_eq!(table.entry(n[2], 0), DofEntry::Fixed);
    }

    #[test]
    fn cycles_and_duplicate_slaves_are_reported() {
        let (map, n) = nodes(2);
        let mut table = DofTable::<NodeId, 3>::for_nodes(map.keys());
        let cycle = vec![
            LinearConstraint {
                slave: (n[0], 0),
                terms: vec![(n[1], 0, 1.0)],
            },
            LinearConstraint {
                slave: (n[1], 0),
                terms: vec![(n[0], 0, 1.0)],
            },
        ];
        assert!(matches!(
            resolve(&cycle, &mut table),
            Err(ResolveError::Cycle(_))
        ));

        let selfref = vec![LinearConstraint {
            slave: (n[0], 0),
            terms: vec![(n[0], 0, 1.0)],
        }];
        assert!(matches!(
            resolve(&selfref, &mut table),
            Err(ResolveError::Cycle(_))
        ));

        let twice = vec![
            LinearConstraint {
                slave: (n[0], 0),
                terms: vec![(n[1], 0, 1.0)],
            },
            LinearConstraint {
                slave: (n[0], 0),
                terms: vec![(n[1], 1, 1.0)],
            },
        ];
        assert_eq!(
            resolve(&twice, &mut table),
            Err(ResolveError::DuplicateSlave((n[0], 0)))
        );
    }
}
