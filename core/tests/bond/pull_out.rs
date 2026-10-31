//! Pull-out of a bar from rigid concrete (every concrete-image node is fixed), loaded at its end. One free slip DOF per
//! rod node: the bar tangent `t`, with the normal motion tied to the image node by a linear constraint (plan §4.1).
//! Units: N, mm, MPa. The reference load is 1 N, so `StepResult::load_factor` is the applied force in N.

use std::f64::consts::PI;

use carapace_core::analysis::{
    Algorithm, Analysis, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
    TangentStrategy,
};
use carapace_core::model::{
    Domain, Element, LimitKind, LimitStateConfig, Material, Node, NodeId, Orientation, Truss,
    ZeroLength,
};

const ES: f64 = 200_000.0;
const EC: f64 = 33_000.0;
const KG: f64 = 0.2;
const RF: f64 = 1e-5;
/// Stop slip `alpha * delta_u_max` (book §3.6.7: alpha = 10; delta_u_max is not given numerically, 1 mm assumed).
const SLIP_LIMIT: f64 = 10.0;

/// A straight bar from the loaded end `x = 0` into the concrete along `angle` (radians from global x).
#[derive(Clone, Copy)]
struct Pull {
    diameter: f64,
    embedment: f64,
    segments: usize,
    fbd: f64,
    angle: f64,
    /// Bond stays elastic (`fy` effectively infinite) when false.
    plastic: bool,
    /// Anchorage spring at the far end, `(Fau, Ka)` (plan §1, book Eq. 3.15-3.16), in parallel with the bond spring there.
    end: Option<(f64, f64)>,
}

impl Pull {
    fn area(&self) -> f64 {
        PI * self.diameter * self.diameter / 4.0
    }
    fn perimeter(&self) -> f64 {
        PI * self.diameter
    }
    fn gb(&self) -> f64 {
        KG * EC / self.diameter
    }
    fn h(&self) -> f64 {
        self.embedment / self.segments as f64
    }
    fn tangent(&self) -> [f64; 2] {
        [self.angle.cos(), self.angle.sin()]
    }
    /// Tributary length of rod node `k` (halves at the ends).
    fn tributary(&self, k: usize) -> f64 {
        if k == 0 || k == self.segments {
            self.h() / 2.0
        } else {
            self.h()
        }
    }
    /// Closed-form pull-out capacity from bond alone.
    fn bond_capacity(&self) -> f64 {
        self.perimeter() * self.embedment * self.fbd
    }
}

struct Model {
    domain: Domain,
    rod: Vec<NodeId>,
    image: Vec<NodeId>,
}

fn build(bar: Pull) -> Model {
    let [tx, ty] = bar.tangent();
    let (nx, ny) = (-ty, tx);
    let mut domain = Domain::new();
    let (mut rod, mut image) = (Vec::new(), Vec::new());
    for k in 0..=bar.segments {
        let (x, y) = (k as f64 * bar.h() * tx, k as f64 * bar.h() * ty);
        image.push(domain.add_node(Node::new([x, y]).fix(0).fix(1).fix(2)));
        rod.push(domain.add_node(Node::new([x, y])));
    }
    for k in 0..=bar.segments {
        // Normal tie: n . (u_R - u_C) = 0, solved for the DOF with the larger normal component.
        let (slave, other, ratio) = if ny.abs() >= nx.abs() {
            (1, 0, nx / ny)
        } else {
            (0, 1, ny / nx)
        };
        domain.add_constraint(
            (rod[k], slave),
            &[
                (image[k], slave, 1.0),
                (rod[k], other, -ratio),
                (image[k], other, ratio),
            ],
        );
        let l = bar.tributary(k);
        let (e0, fy) = (
            bar.gb() * bar.perimeter() * l,
            if bar.plastic {
                bar.fbd * bar.perimeter() * l
            } else {
                1e30
            },
        );
        let spring = ZeroLength::new(image[k], rod[k])
            .with_material(0, Material::steel01(fy, e0, RF, 0.0, 1.0, 0.0, 1.0))
            .with_slip_limit(0, SLIP_LIMIT, LimitKind::BondSlip)
            .with_orientation(Orientation::in_plane(tx, ty).unwrap())
            .unwrap();
        domain.add_element(Element::ZeroLength(spring));
    }
    if let Some((fau, ka)) = bar.end {
        let spring = ZeroLength::new(image[bar.segments], rod[bar.segments])
            .with_material(0, Material::steel01(fau, ka, 1e-2, 0.0, 1.0, 0.0, 1.0))
            .with_slip_limit(0, SLIP_LIMIT, LimitKind::AnchorageSlip)
            .with_orientation(Orientation::in_plane(tx, ty).unwrap())
            .unwrap();
        domain.add_element(Element::ZeroLength(spring));
    }
    for k in 0..bar.segments {
        domain.add_element(Element::Truss(Truss::new(
            rod[k],
            rod[k + 1],
            bar.area(),
            Material::Elastic { e: ES },
        )));
    }
    // The pulling force acts on the loaded end, away from the concrete: along -t. Unit reference load.
    domain.load_node(rod[0], 0, -tx);
    domain.load_node(rod[0], 1, -ty);
    Model { domain, rod, image }
}

fn analysis(domain: Domain, increment: f64, min_fraction: f64) -> Analysis {
    AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControlCutbacks {
            increment,
            min_increment: increment * min_fraction,
        })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::Combined {
            force_tol: 1e-6,
            moment_tol: 1e-6,
            relative_tol: 1e-10,
            displacement_tol: None,
            max_iter: 50,
        })
        .build(domain)
}

/// Loads to `force` in `steps` equal increments; returns the load factor reached (the force, N).
fn load_to(analysis: &mut Analysis, force: f64, steps: usize) -> f64 {
    let mut reached = 0.0;
    for _ in 0..steps {
        reached = analysis
            .step()
            .unwrap_or_else(|e| panic!("stopped at {reached} of {force} N: {e:?}"))
            .load_factor;
    }
    reached
}

/// Slip of rod node `k` along `t` (positive into the concrete): the relative displacement to its fixed image.
fn slip(model: &Model, analysis: &Analysis, bar: &Pull, k: usize) -> f64 {
    let [tx, ty] = bar.tangent();
    let d = analysis.domain().node(model.rod[k]).displacement;
    d[0] * tx + d[1] * ty
}

/// Axial force of segment `e` (tension positive), from the elongation.
fn bar_force(model: &Model, analysis: &Analysis, bar: &Pull, e: usize) -> f64 {
    ES * bar.area() * (slip(model, analysis, bar, e + 1) - slip(model, analysis, bar, e)) / bar.h()
}

fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "{what}: {actual} vs {expected}"
    );
}

/// A1: elastic bond. With the force `F` at the loaded end and a free far end, `N(x) = F sinh(w(l-x))/sinh(wl)` and the
/// slip is `s(x) = F cosh(w(l-x))/(Es As w sinh(wl))`, `w^2 = 4 Gb/(Es D)` (bar pulled out along -t, so the slip
/// along t is negative).
fn elastic_pull_out(angle: f64) {
    let bar = Pull {
        diameter: 16.0,
        embedment: 100.0,
        segments: 40,
        fbd: 3.0,
        angle,
        plastic: false,
        end: None,
    };
    let model = build(bar);
    let mut analysis = analysis(model.domain.clone(), 2_000.0, 1.0 / 16.0);
    let force = load_to(&mut analysis, 10_000.0, 5);
    close(force, 10_000.0, 1e-12, "load reached");

    let w = (4.0 * bar.gb() / (ES * bar.diameter)).sqrt();
    let l = bar.embedment;
    for e in 0..bar.segments {
        let x = (e as f64 + 0.5) * bar.h();
        let expected = force * (w * (l - x)).sinh() / (w * l).sinh();
        let actual = bar_force(&model, &analysis, &bar, e);
        // Compare in the force scale: the far end is near zero, where a relative tolerance is meaningless.
        assert!(
            (actual - expected).abs() <= 1e-2 * force,
            "N at x = {x}: {actual} vs {expected}"
        );
    }
    for k in [0, bar.segments / 2, bar.segments] {
        let x = k as f64 * bar.h();
        let expected = -force * (w * (l - x)).cosh() / (ES * bar.area() * w * (w * l).sinh());
        close(
            slip(&model, &analysis, &bar, k),
            expected,
            1e-2,
            &format!("slip at x = {x}"),
        );
    }
}

#[test]
fn a1_elastic_bond_matches_the_closed_form_along_the_global_axis() {
    elastic_pull_out(0.0);
}

/// The normal tie must not matter: the same bar at 30 degrees gives the same answer.
#[test]
fn a1_elastic_bond_matches_the_closed_form_on_a_skewed_bar() {
    elastic_pull_out(30f64.to_radians());
}

/// A2: rigid-plastic bond at half the capacity. The yielded zone carries `P fbd` per length, so the segment forces
/// equal `F - P fbd x_mid` there (exact for lumped springs at the segment midpoints) and fall to zero beyond the front.
#[test]
fn a2_plastic_bond_carries_p_fbd_per_length_over_the_developed_length() {
    let bar = Pull {
        diameter: 16.0,
        embedment: 500.0,
        segments: 100,
        fbd: 3.0,
        angle: 0.0,
        plastic: true,
        end: None,
    };
    let model = build(bar);
    let target = 0.5 * bar.bond_capacity();
    let mut analysis = analysis(model.domain.clone(), target / 10.0, 1.0 / 64.0);
    close(
        load_to(&mut analysis, target, 10),
        target,
        1e-12,
        "load reached",
    );

    let developed = target / (bar.perimeter() * bar.fbd);
    assert!((developed - 250.0).abs() < 1e-9);
    for e in 0..bar.segments {
        let x = (e as f64 + 0.5) * bar.h();
        let n = bar_force(&model, &analysis, &bar, e);
        if x < developed - 50.0 {
            close(
                n,
                target - bar.perimeter() * bar.fbd * x,
                1e-2,
                &format!("N at x = {x}"),
            );
        }
        if x > developed + 150.0 {
            assert!(
                n.abs() < 1e-2 * target,
                "N beyond the front at x = {x}: {n}"
            );
        }
    }
}

/// Net force the supports (the concrete) hold against the pull, along `t` (N): the sum of the image-node reactions.
fn transferred(model: &Model, analysis: &Analysis, bar: &Pull, load_factor: f64) -> f64 {
    let [tx, ty] = bar.tangent();
    let reactions = analysis.domain().reactions_all(load_factor);
    model
        .image
        .iter()
        .map(|&n| reactions[n][0] * tx + reactions[n][1] * ty)
        .sum::<f64>()
}

/// B4: the bond and anchorage forces carry exactly the applied force into the concrete.
#[test]
fn b4_the_concrete_takes_the_applied_force() {
    for angle in [0.0, 30f64.to_radians()] {
        let bar = Pull {
            diameter: 16.0,
            embedment: 500.0,
            segments: 100,
            fbd: 3.0,
            angle,
            plastic: true,
            end: None,
        };
        let model = build(bar);
        let target = 0.5 * bar.bond_capacity();
        let mut analysis = analysis(model.domain.clone(), target / 10.0, 1.0 / 64.0);
        let force = load_to(&mut analysis, target, 10);
        close(
            transferred(&model, &analysis, &bar, force),
            target,
            1e-6,
            "force into the concrete",
        );
    }
}

/// A3: the capacity is the bond force over the embedment. Below it the slip at the loaded end is the bar stretch over
/// the developed length, `F (f l)/(2 Es As)` with `F = f Fu` (plus the yield slip); just above it the slip runs away,
/// limited only by the 1e-5 plastic slope.
#[test]
fn a3_slip_runs_away_at_the_bond_capacity() {
    let bar = Pull {
        diameter: 16.0,
        embedment: 500.0,
        segments: 100,
        fbd: 3.0,
        angle: 0.0,
        plastic: true,
        end: None,
    };
    close(
        bar.bond_capacity(),
        75_398.2,
        1e-5,
        "book Table 4.6, D16 at 500 mm, fbd 3 MPa: 75.4 kN",
    );
    let capacity = bar.bond_capacity();
    let slip_at = |fraction: f64| {
        let model = build(bar);
        let mut analysis = analysis(
            model.domain.clone(),
            fraction * capacity / 20.0,
            1.0 / 256.0,
        );
        load_to(&mut analysis, fraction * capacity, 20);
        slip(&model, &analysis, &bar, 0).abs()
    };
    let f = 0.98;
    let expected = f * f * capacity * bar.embedment / (2.0 * ES * bar.area());
    let (below, above) = (slip_at(f), slip_at(1.02));
    close(below, expected, 0.05, "slip below the capacity");
    assert!(
        above > 20.0 * below,
        "slip just above the capacity: {above} mm vs {below} mm below"
    );
}

/// E1-E2: an end spring adds `Fau` to what the bar can carry. Without it the slip runs away just above the bond
/// capacity. With it the far end carries force only once the bond is fully developed (F above `P l fbd`), the spring
/// yields at `Fau/Ka`, and beyond that it hardens at `KSN Ka` (so the book stops the run at the slip limit instead).
#[test]
fn e2_the_end_spring_adds_its_force_to_the_capacity() {
    let (diameter, fy, beta) = (16.0, 500.0, 0.3);
    let probe = Pull {
        diameter,
        embedment: 200.0,
        segments: 40,
        fbd: 3.0,
        angle: 0.0,
        plastic: true,
        end: None,
    };
    let area = probe.area();
    let lb_net = diameter * fy / (4.0 * probe.fbd);
    let (fau, ka) = (beta * area * fy, beta * lb_net * KG * EC);
    let bar = Pull {
        end: Some((fau, ka)),
        ..probe
    };
    let (bond, total) = (bar.bond_capacity(), bar.bond_capacity() + fau);
    let slips = |bar: Pull, force: f64| {
        let model = build(bar);
        let mut analysis = analysis(model.domain.clone(), force / 20.0, 1.0 / 256.0);
        load_to(&mut analysis, force, 20);
        (
            slip(&model, &analysis, &bar, 0).abs(),
            slip(&model, &analysis, &bar, bar.segments).abs(),
        )
    };
    assert!(
        slips(probe, 1.02 * bond).0 > 20.0 * slips(probe, 0.98 * bond).0,
        "no end spring: runaway above P l fbd"
    );
    let yield_slip = fau / ka;
    let (_, far_below) = slips(bar, 0.9 * bond);
    // 0.9 of the capacity leaves 20 mm undeveloped, inside the elastic transfer length 1/w = 44 mm, so a small tail remains.
    assert!(
        far_below < 0.1 * yield_slip,
        "the far end is nearly unloaded while the bond is not fully developed: {far_below}"
    );
    // Bond fully developed: the far end carries the rest of the force, F - P l fbd, through the spring.
    let (_, far_at) = slips(bar, 0.99 * total);
    close(
        far_at,
        (0.99 * total - bond) / ka,
        0.03,
        "end-spring slip just below its yield force",
    );
    let (_, far_over) = slips(bar, 1.05 * total);
    close(
        far_over,
        yield_slip + (1.05 * total - bond - fau) / (1e-2 * ka),
        0.05,
        "end-spring slip on its hardening branch",
    );
}

/// U1-U2: bond slip is permanent. Pull the loaded end well past the yield slip and back to zero displacement: a
/// path-independent law would return to zero force and zero slip; the hysteretic one keeps a residual.
#[test]
fn u1_u2_unloading_keeps_the_permanent_slip_and_reverses_at_minus_fbd() {
    let bar = Pull {
        diameter: 16.0,
        embedment: 200.0,
        segments: 40,
        fbd: 3.0,
        angle: 0.0,
        plastic: true,
        end: None,
    };
    let model = build(bar);
    let mut analysis = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::DisplacementControl {
            node: model.rod[0],
            dof: 0,
            increment: 0.0,
        })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::Combined {
            force_tol: 1e-6,
            moment_tol: 1e-6,
            relative_tol: 1e-10,
            displacement_tol: None,
            max_iter: 50,
        })
        .build(model.domain.clone());
    let step = 0.002;
    let slip_yield = bar.fbd / bar.gb();
    let mut peak_force = 0.0f64;
    carapace_core::testkit::ramp_to(&mut analysis, model.rod[0], 0, -0.4, step, |_, f| {
        peak_force = peak_force.max(f.abs())
    });
    assert!(peak_force > 0.9 * bar.bond_capacity(), "peak {peak_force}");
    assert!(0.4 > 20.0 * slip_yield);

    let mut residual = 0.0;
    carapace_core::testkit::ramp_to(&mut analysis, model.rod[0], 0, 0.0, step, |_, f| {
        residual = f
    });
    let mid = bar.segments / 2;
    let permanent = slip(&model, &analysis, &bar, mid).abs();
    assert!(
        permanent > 5.0 * slip_yield,
        "permanent slip at mid-embedment: {permanent} mm"
    );
    assert!(
        residual.abs() > 0.05 * peak_force,
        "residual force at zero end displacement: {residual} N"
    );

    // Reverse: the bond springs yield again at the opposite limit; no spring exceeds its own yield force.
    let mut worst = 0.0f64;
    carapace_core::testkit::ramp_to(&mut analysis, model.rod[0], 0, 0.4, step, |_, f| {
        worst = worst.max(f.abs())
    });
    assert!(worst > 0.9 * bar.bond_capacity(), "reverse peak {worst}");
}

/// A4: the slip limit is a material-limit utilization like concrete crushing: `|s| / (alpha delta_u_max)` of the worst
/// bond spring, reported as `BondSlip` (or `AnchorageSlip` for an end spring), below 1 under the capacity and above it
/// once the bar pulls out.
#[test]
fn a4_the_slip_limit_is_reported_as_a_material_limit() {
    let bar = Pull {
        diameter: 16.0,
        embedment: 500.0,
        segments: 100,
        fbd: 3.0,
        angle: 0.0,
        plastic: true,
        end: None,
    };
    let utilization = |fraction: f64| {
        let model = build(bar);
        let mut analysis = analysis(
            model.domain.clone(),
            fraction * bar.bond_capacity() / 20.0,
            1.0 / 256.0,
        );
        load_to(&mut analysis, fraction * bar.bond_capacity(), 20);
        let (_, utilization, kind) = analysis
            .domain()
            .controlling_limit_utilization(&LimitStateConfig::default())
            .unwrap();
        (utilization, kind, slip(&model, &analysis, &bar, 0).abs())
    };
    let (below, kind, s) = utilization(0.98);
    assert_eq!(kind, LimitKind::BondSlip);
    close(
        below,
        s / SLIP_LIMIT,
        1e-9,
        "utilization = slip / limit at the loaded end",
    );
    assert!(below < 0.1, "{below}");
    let (above, _, _) = utilization(1.02);
    assert!(above > 1.0, "utilization above the capacity: {above}");

    // The end spring reports its own kind.
    let probe = Pull {
        embedment: 200.0,
        segments: 40,
        end: Some((30_000.0, 1.3e6)),
        ..bar
    };
    let model = build(probe);
    let mut analysis = analysis(model.domain.clone(), 60_000.0 / 20.0, 1.0 / 256.0);
    load_to(&mut analysis, 60_000.0, 20);
    let (_, _, kind) = analysis
        .domain()
        .controlling_limit_utilization(&LimitStateConfig::default())
        .unwrap();
    assert!(matches!(
        kind,
        LimitKind::BondSlip | LimitKind::AnchorageSlip
    ));
}
