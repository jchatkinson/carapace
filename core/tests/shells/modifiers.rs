//! ETABS-style stiffness modifiers (`f11 f22 f12 m11 m22 m12 v13 v23`) on `ElasticMembranePlate`.
//!
//! OpenSees has no such modifiers, so nothing here is an OpenSees comparison. The section tests check the
//! tangent against the closed-form `D` the section is documented to build; the element tests use scalings
//! that are exact for the element (a bending-only scaling of `Shell3`, which has no transverse shear, and a
//! uniform scaling of every term of `Shell4`, which equals scaling `E`). Tolerances are stated per test.

use carapace_core::analysis::{
    modal_analysis, Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest,
    Integrator,
};
use carapace_core::model::{
    Domain3, Element3, ElementLoad3, Node3, Node3Id, Shell3, Shell4, ShellModifiers, ShellSection,
    ShellSectionError, ShellVector,
};

const E: f64 = 2.0e5;
const NU: f64 = 0.3;
const H: f64 = 0.2;

fn modified(k: ShellModifiers) -> ShellSection {
    ShellSection::elastic_membrane_plate_modified(E, NU, H, 1.0e-3, k).unwrap()
}

fn close(name: &str, got: f64, want: f64, rel: f64) {
    assert!(
        (got - want).abs() <= rel * want.abs().max(f64::MIN_POSITIVE),
        "{name}: got {got}, want {want}"
    );
}

/// With unequal factors every tangent entry is the unmodified one times its factor, the Poisson coupling
/// terms scale by the geometric mean of the two diagonal factors, bending keeps its negative sign, and the
/// matrix stays symmetric. Tolerance 1e-14 relative.
#[test]
fn the_tangent_scales_each_term_by_its_factor() {
    let k = ShellModifiers {
        f11: 0.4,
        f22: 0.9,
        f12: 0.7,
        m11: 0.25,
        m22: 0.6,
        m12: 0.5,
        v13: 0.3,
        v23: 0.8,
    };
    let base = modified(ShellModifiers::default()).tangent();
    let t = modified(k).tangent();
    let scale = [
        (0, 0, k.f11),
        (1, 1, k.f22),
        (0, 1, (k.f11 * k.f22).sqrt()),
        (1, 0, (k.f11 * k.f22).sqrt()),
        (2, 2, k.f12),
        (3, 3, k.m11),
        (4, 4, k.m22),
        (3, 4, (k.m11 * k.m22).sqrt()),
        (4, 3, (k.m11 * k.m22).sqrt()),
        (5, 5, k.m12),
        (6, 6, k.v13),
        (7, 7, k.v23),
    ];
    for (i, j, factor) in scale {
        close(
            &format!("T[{i}][{j}]"),
            t[(i, j)],
            base[(i, j)] * factor,
            1e-14,
        );
    }
    assert!(t[(3, 3)] < 0.0 && t[(4, 4)] < 0.0 && t[(5, 5)] < 0.0);
    assert!((t - t.transpose()).abs().max() == 0.0);
}

/// A pure `kappa_x` curvature gives `Mx` scaled by `m11` and `My` by the geometric mean `sqrt(m11 m22)`;
/// a mixed strain field matches the full formula, so `Mx` is not simply `m11` times the unmodified value.
#[test]
fn resultants_follow_the_modified_bending_block() {
    let k = ShellModifiers {
        m11: 0.25,
        m22: 0.64,
        ..ShellModifiers::default()
    };
    let d = E * H.powi(3) / 12.0 / (1.0 - NU * NU);
    let section = modified(k);
    let mut strain = ShellVector::zeros();
    strain[3] = 1e-3;
    let (r, _) = section.trial_resultant(&strain);
    close("Mx", r[3], -d * k.m11 * 1e-3, 1e-14);
    close("My", r[4], -NU * d * (k.m11 * k.m22).sqrt() * 1e-3, 1e-14);
    strain[4] = 2e-3;
    let (r, _) = section.trial_resultant(&strain);
    let mx = -d * (k.m11 * 1e-3 + NU * (k.m11 * k.m22).sqrt() * 2e-3);
    close("Mx mixed", r[3], mx, 1e-14);
}

/// All factors 1 leave the tangent and the drilling stiffness identical to the unmodified section.
#[test]
fn unit_modifiers_change_nothing() {
    let plain = ShellSection::elastic_membrane_plate(E, NU, H, 1.0e-3).unwrap();
    let unit = modified(ShellModifiers::default());
    assert_eq!(plain.tangent(), unit.tangent());
    assert_eq!(plain.drilling_stiffness(), unit.drilling_stiffness());
}

/// The drilling stiffness is the smallest eigenvalue of the modified membrane block: a uniform factor scales
/// it, and extreme ratios neither return zero nor underflow. Tolerance 1e-12 relative.
#[test]
fn the_drilling_stiffness_follows_the_modified_membrane_block() {
    let m = E / (1.0 - NU * NU) * H;
    let base = modified(ShellModifiers::default()).drilling_stiffness();
    let c = 0.37;
    let uniform = modified(ShellModifiers {
        f11: c,
        f22: c,
        f12: c,
        ..ShellModifiers::default()
    });
    close("uniform", uniform.drilling_stiffness(), c * base, 1e-12);
    // f11 = 1e-20, f22 = 1: the smaller eigenvalue is about m f11 (1 - nu^2), far below the shear term.
    let thin = modified(ShellModifiers {
        f11: 1e-20,
        ..ShellModifiers::default()
    });
    close(
        "1e-20",
        thin.drilling_stiffness(),
        m * 1e-20 * (1.0 - NU * NU),
        1e-9,
    );
    let tiny = modified(ShellModifiers {
        f11: 1e-200,
        f22: 1e-200,
        f12: 1e-200,
        ..ShellModifiers::default()
    });
    let kd = tiny.drilling_stiffness();
    assert!(kd.is_finite() && kd > 0.0, "drilling stiffness {kd}");
    close("1e-200", kd, 1e-200 * base, 1e-9);
}

/// Zero, negative and non-finite factors are rejected.
#[test]
fn invalid_modifiers_are_rejected() {
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        for field in 0..8 {
            let mut k = ShellModifiers::default();
            *[
                &mut k.f11, &mut k.f22, &mut k.f12, &mut k.m11, &mut k.m22, &mut k.m12, &mut k.v13,
                &mut k.v23,
            ][field] = bad;
            assert_eq!(
                ShellSection::elastic_membrane_plate_modified(E, NU, H, 1e-3, k),
                Err(ShellSectionError::InvalidModifier),
                "factor {field} = {bad}"
            );
        }
    }
}

/// A simply supported `a x a` plate of `n x n` cells, each a `Shell4` or two `Shell3`s, under pressure.
/// Returns the domain and the centre node.
fn plate(triangles: bool, section: &ShellSection, pressure: f64, a: f64) -> (Domain3, Node3Id) {
    let n = 8;
    let mut domain = Domain3::new();
    let mut at = vec![];
    for j in 0..=n {
        for i in 0..=n {
            let (x, y) = (a * i as f64 / n as f64, a * j as f64 / n as f64);
            let mut node = Node3::new([x, y, 0.0]);
            if i == 0 || j == 0 || i == n || j == n {
                node = node.fix(2);
            }
            if i == 0 && j == 0 {
                node = node.fix(0).fix(1);
            }
            if i == n && j == 0 {
                node = node.fix(1);
            }
            at.push(domain.add_node(node));
        }
    }
    let pattern = domain.default_pattern();
    for j in 0..n {
        for i in 0..n {
            let q = [
                at[j * (n + 1) + i],
                at[j * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i + 1],
                at[(j + 1) * (n + 1) + i],
            ];
            let cells: Vec<Element3> = if triangles {
                vec![
                    Element3::Shell3(Shell3::new([q[0], q[1], q[2]], section.clone())),
                    Element3::Shell3(Shell3::new([q[0], q[2], q[3]], section.clone())),
                ]
            } else {
                vec![Element3::Shell4(Shell4::new(q, section.clone()))]
            };
            for cell in cells {
                let id = domain.add_element(cell);
                if pressure != 0.0 {
                    domain.add_element_load(pattern, id, ElementLoad3::pressure(pressure));
                }
            }
        }
    }
    let centre = at[(n / 2) * (n + 1) + n / 2];
    (domain, centre)
}

fn solve(domain: Domain3) -> Analysis<3, 6, Node3Id, Element3> {
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Plain)
        .integrator(Integrator::LoadControl { increment: 1.0 })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-9,
            max_iter: 5,
        })
        .build(domain);
    analysis.step().expect("linear solve");
    analysis
}

fn centre_deflection(triangles: bool, section: &ShellSection, a: f64) -> f64 {
    let (domain, centre) = plate(triangles, section, 0.01, a);
    solve(domain).domain().node(centre).displacement[2]
}

fn fundamental_frequency(triangles: bool, section: &ShellSection) -> f64 {
    let (mut domain, _) = plate(triangles, section, 0.0, 10.0);
    modal_analysis(&mut domain, 1).expect("modal analysis")[0].frequency
}

/// `Shell3` has no transverse shear and its bending decouples from the membrane, so scaling all three bending
/// factors by `k` scales a plate's deflection by exactly `1 / k`. Tolerance 1e-9 relative.
#[test]
fn shell3_bending_modifiers_scale_the_deflection_exactly() {
    let k = 0.35;
    let base = centre_deflection(true, &modified(ShellModifiers::default()), 10.0);
    let soft = centre_deflection(
        true,
        &modified(ShellModifiers {
            m11: k,
            m22: k,
            m12: k,
            ..ShellModifiers::default()
        }),
        10.0,
    );
    close("deflection", soft, base / k, 1e-9);
}

/// With every factor `c` the `Shell4` plate equals the same plate with `E * c`: deflection `1 / c`, and the
/// fundamental frequency scales by `sqrt(c)` because the mass is unmodified. Tolerance 1e-9 relative.
#[test]
fn shell4_uniform_modifiers_equal_scaling_the_modulus() {
    let c = 0.4;
    let uniform = modified(ShellModifiers {
        f11: c,
        f22: c,
        f12: c,
        m11: c,
        m22: c,
        m12: c,
        v13: c,
        v23: c,
    });
    let scaled_e = ShellSection::elastic_membrane_plate(E * c, NU, H, 1.0e-3).unwrap();
    close(
        "deflection",
        centre_deflection(false, &uniform, 10.0),
        centre_deflection(false, &scaled_e, 10.0),
        1e-9,
    );
    let base = fundamental_frequency(false, &modified(ShellModifiers::default()));
    close(
        "frequency",
        fundamental_frequency(false, &uniform),
        base * c.sqrt(),
        1e-9,
    );
}

/// `Shell3` modal check: bending factor `k` scales the fundamental plate frequency by `sqrt(k)`.
/// Tolerance 1e-9 relative.
#[test]
fn shell3_bending_modifiers_scale_the_frequency() {
    let k = 0.5;
    let base = fundamental_frequency(true, &modified(ShellModifiers::default()));
    let soft = fundamental_frequency(
        true,
        &modified(ShellModifiers {
            m11: k,
            m22: k,
            m12: k,
            ..ShellModifiers::default()
        }),
    );
    close("frequency", soft, base * k.sqrt(), 1e-9);
}

/// A single `Shell4` strip under an end load in its plane: with free sides the axial stiffness is
/// `E h f11 / L * H` whatever `f22`, so `f11 = k` scales the extension by `1 / k` and `f22` leaves it alone.
/// Tolerance 1e-9 relative.
#[test]
fn the_membrane_factor_f11_scales_the_axial_stiffness() {
    let extension = |k: ShellModifiers| {
        let (l, w) = (3.0, 1.0);
        let mut domain = Domain3::new();
        let coords = [[0.0, 0.0], [l, 0.0], [l, w], [0.0, w]];
        let ids: Vec<Node3Id> = coords
            .iter()
            .enumerate()
            .map(|(i, &[x, y])| {
                // Keep bending and the out-of-plane motion out of the problem.
                let mut node = Node3::new([x, y, 0.0]).fix(2).fix(3).fix(4);
                if x == 0.0 {
                    node = node.fix(0);
                }
                if i == 0 {
                    node = node.fix(1);
                }
                domain.add_node(node)
            })
            .collect();
        let q = [ids[0], ids[1], ids[2], ids[3]];
        domain.add_element(Element3::Shell4(Shell4::new(q, modified(k))));
        domain.load_node(ids[1], 0, 0.5);
        domain.load_node(ids[2], 0, 0.5);
        solve(domain).domain().node(ids[1]).displacement[0]
    };
    let base = extension(ShellModifiers::default());
    close("base", base, 3.0 / (E * H * 1.0), 1e-9);
    let stiff = extension(ShellModifiers {
        f11: 2.5,
        ..ShellModifiers::default()
    });
    close("f11", stiff, base / 2.5, 1e-9);
    let other = extension(ShellModifiers {
        f22: 0.3,
        ..ShellModifiers::default()
    });
    close("f22", other, base, 1e-9);
}

/// Transverse shear factors soften a thick `Shell4` plate, and do nothing to `Shell3`, whose shear rows are
/// zero. Tolerance for the `Shell3` equality 1e-12 relative.
#[test]
fn shear_modifiers_act_on_shell4_and_not_on_shell3() {
    let thick = |v: f64| {
        ShellSection::elastic_membrane_plate_modified(
            E,
            NU,
            2.0,
            1.0e-3,
            ShellModifiers {
                v13: v,
                v23: v,
                ..ShellModifiers::default()
            },
        )
        .unwrap()
    };
    let rigid = centre_deflection(false, &thick(1.0), 10.0);
    let soft = centre_deflection(false, &thick(0.05), 10.0);
    assert!(soft > 1.01 * rigid, "soft {soft} vs {rigid}");
    close(
        "shell3",
        centre_deflection(true, &thick(0.05), 10.0),
        centre_deflection(true, &thick(1.0), 10.0),
        1e-12,
    );
}
