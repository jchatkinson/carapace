use slotmap::{Key, SecondaryMap};

/// What one `(node, dof)` slot contributes to the free-DOF system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DofEntry {
    /// No equation: a fixed DOF, a DOF identity-tied to a fixed one, or an
    /// affine tie whose every term is fixed.
    Fixed,
    /// A free unknown (or a DOF identity-tied to one: tied DOFs share the
    /// retained DOF's equation number).
    Free(usize),
    /// A DOF that is an affine combination of free unknowns:
    /// `value = sum(coeff * q[eq])` over `DofTable::terms(start, len)`.
    Affine { start: u32, len: u32 },
}

/// The `Domain`'s record of how every node DOF maps onto the free-DOF
/// system, produced by `Domain::number_dofs`. Dense (indexed by a per-node
/// ordinal, then DOF slot) rather than a `HashMap`, because assembly and
/// state scatter/gather read it in their hot loops. Replaces the
/// per-node `equation` array and the affine-term `HashMap` the domain
/// previously kept.
///
/// A node that was added after the last numbering has no ordinal; every
/// lookup for it reports `DofEntry::Fixed`, matching the old behavior of a
/// freshly added node having no equation numbers until the next
/// `number_dofs`.
#[derive(Debug, Clone)]
pub(crate) struct DofTable<NId: Key, const NDOF: usize> {
    ordinal: SecondaryMap<NId, usize>,
    entries: Vec<DofEntry>,
    terms: Vec<(usize, f64)>,
}

impl<NId: Key, const NDOF: usize> Default for DofTable<NId, NDOF> {
    fn default() -> Self {
        DofTable {
            ordinal: SecondaryMap::new(),
            entries: Vec::new(),
            terms: Vec::new(),
        }
    }
}

impl<NId: Key, const NDOF: usize> DofTable<NId, NDOF> {
    /// A table covering `nodes` (in iteration order), every DOF `Fixed`.
    pub(crate) fn for_nodes(nodes: impl Iterator<Item = NId>) -> Self {
        let mut ordinal = SecondaryMap::new();
        let mut count = 0;
        for id in nodes {
            ordinal.insert(id, count);
            count += 1;
        }
        DofTable {
            ordinal,
            entries: vec![DofEntry::Fixed; count * NDOF],
            terms: Vec::new(),
        }
    }

    pub(crate) fn entry(&self, node: NId, dof: usize) -> DofEntry {
        assert!(dof < NDOF, "DOF index {dof} out of range for a {NDOF}-DOF node");
        match self.ordinal.get(node) {
            Some(&ordinal) => self.entries[ordinal * NDOF + dof],
            None => DofEntry::Fixed,
        }
    }

    pub(crate) fn set(&mut self, node: NId, dof: usize, entry: DofEntry) {
        assert!(dof < NDOF, "DOF index {dof} out of range for a {NDOF}-DOF node");
        let ordinal = self.ordinal[node];
        self.entries[ordinal * NDOF + dof] = entry;
    }

    /// Stores an affine term list and returns the entry that refers to it
    /// (`Fixed` if every term was dropped, since such a DOF has no equation
    /// and contributes nothing).
    pub(crate) fn push_terms(&mut self, terms: &[(usize, f64)]) -> DofEntry {
        if terms.is_empty() {
            return DofEntry::Fixed;
        }
        let start = self.terms.len() as u32;
        self.terms.extend_from_slice(terms);
        DofEntry::Affine {
            start,
            len: terms.len() as u32,
        }
    }

    pub(crate) fn terms(&self, start: u32, len: u32) -> &[(usize, f64)] {
        &self.terms[start as usize..(start + len) as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::NodeId;
    use slotmap::SlotMap;

    #[test]
    fn unknown_nodes_and_unset_slots_are_fixed() {
        let mut nodes: SlotMap<NodeId, ()> = SlotMap::with_key();
        let a = nodes.insert(());
        let mut table = DofTable::<NodeId, 3>::for_nodes(nodes.keys());
        assert_eq!(table.entry(a, 1), DofEntry::Fixed);

        table.set(a, 1, DofEntry::Free(7));
        assert_eq!(table.entry(a, 1), DofEntry::Free(7));
        assert_eq!(table.entry(a, 0), DofEntry::Fixed);

        // Added after the table was built: no ordinal, so reads as fixed.
        let b = nodes.insert(());
        assert_eq!(table.entry(b, 0), DofEntry::Fixed);
    }

    #[test]
    fn affine_terms_round_trip_and_empty_terms_collapse_to_fixed() {
        let mut nodes: SlotMap<NodeId, ()> = SlotMap::with_key();
        let a = nodes.insert(());
        let mut table = DofTable::<NodeId, 3>::for_nodes(nodes.keys());
        let entry = table.push_terms(&[(2, 1.0), (5, -0.5)]);
        table.set(a, 0, entry);
        let DofEntry::Affine { start, len } = table.entry(a, 0) else {
            panic!("expected affine entry");
        };
        assert_eq!(table.terms(start, len), &[(2, 1.0), (5, -0.5)]);
        assert_eq!(table.push_terms(&[]), DofEntry::Fixed);
    }
}
