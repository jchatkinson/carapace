//! Material probe cases (unit zero-length spring, prescribed strain), as used by the PySees material preview.

use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::material_probe::{MaterialProbe, MaterialProbeConfig, MaterialProbeError};

fn probe(material: MaterialSpec) -> MaterialProbe {
    MaterialProbe::new(&MaterialProbeConfig {
        material,
        arena: vec![],
        initial_strain: None,
    })
    .unwrap()
}

fn run(probe: &mut MaterialProbe, targets: &[f64]) -> Vec<f64> {
    targets
        .iter()
        .map(|&t| {
            let r = probe.apply_strain(t).unwrap();
            assert_eq!(r.strain, t);
            r.stress
        })
        .collect()
}

fn close(got: &[f64], want: &[f64]) {
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() < 1e-6, "{got:?} vs {want:?}");
    }
}

/// An elastic probe follows a prescribed strain history through zero, positive and negative
/// strains: stress equals E times strain.
#[test]
fn elastic() {
    let mut p = probe(MaterialSpec::Elastic { e: 200_000.0 });
    close(
        &run(&mut p, &[0.0, 0.001, -0.001, 0.0]),
        &[0.0, 200.0, -200.0, 0.0],
    );
}

/// A Concrete01 probe has no tension capacity (zero stress at positive strain), goes into
/// compression with softening past the peak, and unloads to zero stress at zero strain.
#[test]
fn concrete01_tension_then_compression() {
    let mut p = probe(MaterialSpec::Concrete01 {
        fpc: -30.0,
        epsc0: -0.002,
        fpcu: -6.0,
        epscu: -0.006,
    });
    let s = run(&mut p, &[0.0, 0.001, -0.0005, -0.002, 0.0]);
    assert_eq!(s[0], 0.0);
    assert_eq!(s[1], 0.0, "no tension capacity");
    assert!(s[2] < 0.0 && s[3] < s[2], "{s:?}");
    assert!(
        s[4].abs() < 1e-6,
        "unloads to zero stress at zero strain: {s:?}"
    );
}

/// An elastic-perfectly-plastic probe keeps its committed plastic strain across targets: strains of
/// +0.002, -0.002 and -0.0005 give stresses of +100, -100 and +50 (elastic unloading from the
/// permanent set).
#[test]
fn elastic_pp_keeps_committed_plastic_history() {
    let mut p = probe(MaterialSpec::ElasticPp {
        e: 100_000.0,
        eyp: 0.001,
    });
    close(
        &run(&mut p, &[0.002, -0.002, -0.0005]),
        &[100.0, -100.0, 50.0],
    );
}

/// Applying the same target strain repeatedly does not advance material state: the
/// elastic-perfectly-plastic probe returns the same stress each time.
#[test]
fn repeated_targets_do_not_advance_state() {
    let mut p = probe(MaterialSpec::ElasticPp {
        e: 100_000.0,
        eyp: 0.001,
    });
    close(
        &run(&mut p, &[0.0, 0.0, 0.0005, 0.0005]),
        &[0.0, 0.0, 50.0, 50.0],
    );
}

/// `reset` discards the history, and the configured initial strain is applied: with an initial
/// strain of 0.5, a target of 1.0 and later 0.5 give stresses 10 and 5 for E = 10.
#[test]
fn reset_discards_history_and_initial_strain_is_applied() {
    let mut p = MaterialProbe::new(&MaterialProbeConfig {
        material: MaterialSpec::Elastic { e: 10.0 },
        arena: vec![],
        initial_strain: Some(0.5),
    })
    .unwrap();
    close(&run(&mut p, &[1.0]), &[10.0]);
    p.reset().unwrap();
    close(&run(&mut p, &[0.5]), &[5.0]);
}

/// Probe errors are structured: an unsupported material (`Gap`) is `UnsupportedMaterial`, a NaN
/// target is `invalidTarget`, and a rejected target leaves the committed state unchanged.
#[test]
fn errors_are_structured() {
    let unsupported = MaterialProbe::new(&MaterialProbeConfig {
        material: MaterialSpec::Gap { e: 1.0, gap: 0.1 },
        arena: vec![],
        initial_strain: None,
    });
    assert_eq!(
        unsupported.err(),
        Some(MaterialProbeError::UnsupportedMaterial { material: "gap" })
    );
    let mut p = probe(MaterialSpec::Elastic { e: 1.0 });
    p.apply_strain(0.25).unwrap();
    assert_eq!(
        p.apply_strain(f64::NAN).unwrap_err().to_string_kind(),
        "invalidTarget"
    );
    // Rejected target leaves committed state intact.
    assert_eq!(p.apply_strain(0.25).unwrap().stress, 0.25);
}

trait Kind {
    fn to_string_kind(&self) -> &'static str;
}
impl Kind for MaterialProbeError {
    fn to_string_kind(&self) -> &'static str {
        match self {
            MaterialProbeError::InvalidTarget { .. } => "invalidTarget",
            _ => "other",
        }
    }
}

/// Composite materials refer to arena entries by index; the probed material may be any of them, and a
/// `Gap` anywhere in the arena is still rejected.
#[test]
fn composites_resolve_through_the_arena() {
    let config = |material, arena| MaterialProbeConfig {
        material,
        arena,
        initial_strain: None,
    };
    let mut p = MaterialProbe::new(&config(
        MaterialSpec::Parallel {
            children: vec![(0, 1.0), (1, 2.0)],
        },
        vec![
            MaterialSpec::Elastic { e: 10.0 },
            MaterialSpec::Elastic { e: 5.0 },
        ],
    ))
    .unwrap();
    close(&run(&mut p, &[0.1]), &[2.0]);
    let mut p = MaterialProbe::new(&config(
        MaterialSpec::Series {
            children: vec![0, 1],
        },
        vec![
            MaterialSpec::Elastic { e: 10.0 },
            MaterialSpec::Elastic { e: 10.0 },
        ],
    ))
    .unwrap();
    close(&run(&mut p, &[1.0]), &[5.0]);
    let gapped = MaterialProbe::new(&config(
        MaterialSpec::MinMax {
            inner: 0,
            min_strain: -1.0,
            max_strain: 1.0,
        },
        vec![MaterialSpec::Gap { e: 1.0, gap: 0.1 }],
    ));
    assert_eq!(
        gapped.err(),
        Some(MaterialProbeError::UnsupportedMaterial { material: "gap" })
    );
}

/// Steel02 and Concrete02 are probe-able now.
#[test]
fn steel02_and_concrete02_probe() {
    let mut p = probe(MaterialSpec::Steel02 {
        fy: 50.0,
        e0: 1000.0,
        b: 0.01,
        r0: 15.0,
        cr1: 0.925,
        cr2: 0.15,
        a1: 0.0,
        a2: 1.0,
        a3: 0.0,
        a4: 1.0,
    });
    assert!(run(&mut p, &[0.01, 0.0, -0.01])
        .iter()
        .all(|s| s.abs() <= 60.0));
    let mut p = probe(MaterialSpec::Concrete02 {
        fc: -30.0,
        epsc0: -0.002,
        fcu: -20.0,
        epscu: -0.006,
        rat: 0.1,
        ft: 3.0,
        ets: 200.0,
    });
    let r = run(&mut p, &[0.0005, -0.002]);
    assert!(r[0] > 0.0 && r[1] < -20.0);
}
