/// Static/pseudo-static integration strategy. Closed enum (§3.1).
/// `DisplacementControl` lands at M4.
#[derive(Debug, Clone, Copy)]
pub enum Integrator {
    LoadControl { increment: f64 },
}

impl Integrator {
    pub(crate) fn next_load_factor(&self, current: f64) -> f64 {
        match self {
            Integrator::LoadControl { increment } => current + increment,
        }
    }
}
