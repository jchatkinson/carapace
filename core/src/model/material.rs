/// Uniaxial material catalog. Closed enum, not a trait object (§3.1 of the
/// implementation plan) — the set is fixed at compile time. Leaf variants
/// only for M1; recursive composites (Parallel/Series/MinMax) land at M7.
#[derive(Debug, Clone, Copy)]
pub enum Material {
    Elastic { e: f64 },
}

impl Material {
    /// Given the current strain, return `(stress, tangent)`. Leaf materials
    /// are stateless for now (linear elastic has no history to track);
    /// stateful leaves (Steel01, Concrete01, ...) will need `&mut self` and
    /// internal history variables when they land at M7.
    pub fn stress_tangent(&self, strain: f64) -> (f64, f64) {
        match self {
            Material::Elastic { e } => (e * strain, *e),
        }
    }
}
