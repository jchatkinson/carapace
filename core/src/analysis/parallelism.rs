//! Process-wide thread count for the sparse factorizations (native only, `parallel` feature).
//!
//! `faer` reads one global setting, so this is global too: it applies to every analysis in the process.
//! Callers that run many independent analyses concurrently (`carapace_input::run_many`) set `Serial` for
//! the duration so the two levels of parallelism do not oversubscribe the machine. Results do not depend
//! on the setting; `tests/analysis/parallel.rs` checks that bit for bit.

/// How many threads the linear solver may use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Parallelism {
    Serial,
    /// This many threads; `0` means every available core.
    Threads(usize),
}

/// Sets the solver's thread count. Without the `parallel` feature only `Serial` is accepted and anything
/// else returns `false` (the setting is unchanged), so callers can report that the build is serial.
pub fn set_parallelism(parallelism: Parallelism) -> bool {
    #[cfg(feature = "parallel")]
    {
        *POOL.lock().expect("pool lock poisoned") = match parallelism {
            Parallelism::Serial => None,
            Parallelism::Threads(n) => rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .ok()
                .map(std::sync::Arc::new),
        };
        faer::set_global_parallelism(match parallelism {
            Parallelism::Serial => faer::Par::Seq,
            Parallelism::Threads(n) => faer::Par::rayon(n),
        });
        true
    }
    #[cfg(not(feature = "parallel"))]
    {
        parallelism == Parallelism::Serial
    }
}

/// The current setting.
pub fn parallelism() -> Parallelism {
    #[cfg(feature = "parallel")]
    {
        match faer::get_global_parallelism() {
            faer::Par::Seq => Parallelism::Serial,
            faer::Par::Rayon(n) => Parallelism::Threads(n.get()),
        }
    }
    #[cfg(not(feature = "parallel"))]
    {
        Parallelism::Serial
    }
}

/// The pool the element loop runs on, `None` when serial.
#[cfg(feature = "parallel")]
static POOL: std::sync::Mutex<Option<std::sync::Arc<rayon::ThreadPool>>> =
    std::sync::Mutex::new(None);

#[cfg(feature = "parallel")]
pub(crate) fn pool() -> Option<std::sync::Arc<rayon::ThreadPool>> {
    POOL.lock().expect("pool lock poisoned").clone()
}
