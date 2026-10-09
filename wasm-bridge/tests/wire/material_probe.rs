//! Material probe cases (unit zero-length spring, prescribed strain), as used by the PySees material preview.

use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::material_probe::{MaterialProbe, MaterialProbeConfig, MaterialProbeError};

fn probe(material: MaterialSpec) -> MaterialProbe {
    MaterialProbe::new(&MaterialProbeConfig {
        material,
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

#[test]
fn elastic() {
    let mut p = probe(MaterialSpec::Elastic { e: 200_000.0 });
    close(
        &run(&mut p, &[0.0, 0.001, -0.001, 0.0]),
        &[0.0, 200.0, -200.0, 0.0],
    );
}

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

#[test]
fn reset_discards_history_and_initial_strain_is_applied() {
    let mut p = MaterialProbe::new(&MaterialProbeConfig {
        material: MaterialSpec::Elastic { e: 10.0 },
        initial_strain: Some(0.5),
    })
    .unwrap();
    close(&run(&mut p, &[1.0]), &[10.0]);
    p.reset().unwrap();
    close(&run(&mut p, &[0.5]), &[5.0]);
}

#[test]
fn errors_are_structured() {
    let unsupported = MaterialProbe::new(&MaterialProbeConfig {
        material: MaterialSpec::Gap { e: 1.0, gap: 0.1 },
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
