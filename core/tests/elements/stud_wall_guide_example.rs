//! The usage example from `docs/stud-wall.md`, compiled and run so the guide cannot drift from the API.
//! It checks statics only (base shear and vertical reaction equal the applied loads, pier shears sum to the
//! applied shear, M/jd equals the base overturning over `jd`), to 1e-8 relative.

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator, TangentStrategy,
};
use carapace_core::model::{
    ChordDistribution, Domain, Element, ElementLoad, LoadSeries, Material, Node, NodeId, StudWall,
    StudWallLayout, WallLoad, WallOpening, WallShearLaw, WallSheathing, WallStud, ZeroLength,
    WALL_SLOTS,
};

/// Guide example: 3.6 m wall on a foundation, window, O86 sheathing law, hold-down spring, gravity plus
/// wind in one Newton run. Statics: base shear = applied shear and pier `chords` = `∓moment_base/jd`, 1e-8.
#[test]
fn guide_example_runs() {
    // --- guide example start ---
    let (l, h) = (3600.0, 2440.0); // mm
    let ea_ply = 9500.0 * 38.0 * 140.0; // N, one 38x140 SPF ply
    let mut studs: Vec<WallStud> = (0..=9)
        .map(|i| WallStud {
            x: i as f64 * 400.0,
            ea: ea_ply,
        })
        .collect();
    studs[0].ea = 3.0 * ea_ply; // 3-ply end posts
    studs[9].ea = 3.0 * ea_ply;
    // The window's jacks occupy x0 = 1600 and x1 = 2400, so no stud row may sit there. King studs beside the
    // jacks end each pier at the opening; without them the piers would stop at the last regular stud.
    studs.retain(|s| s.x != 1600.0 && s.x != 2400.0);
    studs.push(WallStud {
        x: 1562.0,
        ea: ea_ply,
    });
    studs.push(WallStud {
        x: 2438.0,
        ea: ea_ply,
    });
    studs.sort_by(|a, b| a.x.total_cmp(&b.x));
    let layout = StudWallLayout {
        length: l,
        height: h,
        studs,
        top_ea: 2.0 * ea_ply, // double top plate
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: ea_ply,
        bottom_ei: 9500.0 * 140.0 * 38f64.powi(3) / 12.0,
        openings: vec![WallOpening {
            x0: 1600.0,
            x1: 2400.0,
            head_height: 2000.0,
            header_ei: 9500.0 * 2.0 * 38.0 * 235f64.powi(3) / 12.0,
            jack_ea: ea_ply,
            door: false,
        }],
        sheathing: Some(WallSheathing {
            law: WallShearLaw::O86 {
                bv: 11_000.0,
                nail_spacing: 100.0,
                nail_diameter: 3.25,
            },
        }),
        foundation: true,
    };

    let mut d: Domain<2, 3, NodeId, Element> = Domain::new();
    let ground = d.add_node(Node::new([0.0, 0.0]).fix(0).fix(1).fix(2));
    let c0 = d.add_node(Node::new([0.0, 0.0]).fix(0)); // on the hold-down
    let c1 = d.add_node(Node::new([l, 0.0]).fix(0).fix(1)); // bears on the foundation
    let c2 = d.add_node(Node::new([l, h]));
    let c3 = d.add_node(Node::new([0.0, h]));
    d.add_element(Element::ZeroLength(
        ZeroLength::new(ground, c0).with_material(1, Material::Elastic { e: 20_000.0 }),
    ));

    // Slot 0: roof dead load. Slot 1: wind shear along the top plate.
    let mut loads: [Vec<WallLoad>; WALL_SLOTS] = std::array::from_fn(|_| Vec::new());
    loads[0].push(WallLoad::TopLine {
        x0: 0.0,
        x1: l,
        wx0: 0.0,
        wx1: 0.0,
        wy0: -5.0,
        wy1: -5.0,
    });
    loads[1].push(WallLoad::TopLine {
        x0: 0.0,
        x1: l,
        wx0: 4.0,
        wx1: 4.0,
        wy0: 0.0,
        wy1: 0.0,
    });
    let wall = StudWall::new([c0, c1, c2, c3], vec![], layout, loads).expect("valid layout");
    let id = d.add_element(Element::StudWall(wall));

    for (slot, factor) in [(0, 1.0), (1, 1.0)] {
        let pattern = d.add_load_pattern_scaled(LoadSeries::Linear { slope: 1.0 }, factor);
        let mut wall = [0.0; WALL_SLOTS];
        wall[slot] = 1.0;
        d.add_element_load(
            pattern,
            id,
            ElementLoad {
                wall,
                ..Default::default()
            },
        );
    }

    // Nonlinear sheathing law: Newton, in load steps (load cases do not superpose).
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 0.1 })
        .algorithm(Algorithm::Newton {
            tangent: TangentStrategy::Current,
            line_search: None,
        })
        .test(ConvergenceTest::NormDispIncr {
            tol: 1e-10,
            max_iter: 25,
        })
        .build(d);
    for _ in 0..10 {
        a.step().expect("converged step");
    }
    let r = a.domain().element_stud_wall_response(id, 1.0).unwrap();
    for pier in &r.piers {
        // q (N/mm), M/jd chord forces (N, tension positive) for the hold-down and post design.
        println!(
            "pier {}..{}: q = {:.2}, chords = {:?}",
            pier.xa, pier.xb, pier.q, pier.chords
        );
    }
    // --- guide example end ---

    let applied = 4.0 * l;
    assert!(
        (r.base[0].abs() - applied).abs() <= 1e-8 * applied,
        "base shear {}",
        r.base[0]
    );
    // Gravity: the supports push up on the wall (positive `fy`), and they carry the whole roof load.
    let gravity = 5.0 * l;
    assert!(
        (r.base[1] - gravity).abs() <= 1e-8 * gravity,
        "base Fy {}",
        r.base[1]
    );
    let shear: f64 = r.piers.iter().map(|p| p.q * (p.xb - p.xa)).sum();
    assert!(
        (shear - applied).abs() <= 1e-8 * applied,
        "pier shear {shear}"
    );
    for p in &r.piers {
        let t = -p.moment_base / p.jd;
        assert!((p.chords[0] - t).abs() <= 1e-8 * t.abs().max(1.));
    }
}
