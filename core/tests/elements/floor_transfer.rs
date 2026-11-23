//! One-way floor reactions against hand statics and the Euler–Bernoulli three-moment solution (1e-12 relative).

use carapace_core::analysis::{
    Algorithm, AnalysisBuilder, ConstraintHandler, ConvergenceTest, Integrator,
};
use carapace_core::model::{
    transfer_floor, ChordDistribution, Domain, Element, ElementLoad, FloorBounds, FloorContinuity,
    FloorLoad, FloorSupport, Node, NodeId, StudWall, StudWallLayout, WallLoad, WallOpening,
    WallStud, WALL_SLOTS,
};

fn bounds() -> FloorBounds {
    FloorBounds {
        s0: 0.,
        s1: 10.,
        t0: 0.,
        t1: 4.,
    }
}
fn supports() -> [FloorSupport; 3] {
    [
        FloorSupport {
            s: 0.,
            t0: 0.,
            t1: 4.,
        },
        FloorSupport {
            s: 4.,
            t0: 0.,
            t1: 4.,
        },
        FloorSupport {
            s: 10.,
            t0: 0.,
            t1: 4.,
        },
    ]
}
fn sum(r: &carapace_core::model::FloorReactions, slot: usize) -> f64 {
    r.per_support
        .iter()
        .map(|s| {
            s.lines
                .iter()
                .filter(|x| x.slot == slot)
                .map(|x| x.w * (x.t1 - x.t0))
                .sum::<f64>()
                + s.points
                    .iter()
                    .filter(|x| x.slot == slot)
                    .map(|x| x.p)
                    .sum::<f64>()
        })
        .sum()
}

fn equilibrium(
    b: FloorBounds,
    supports: &[FloorSupport],
    loads: &[FloorLoad],
    r: &carapace_core::model::FloorReactions,
) {
    for slot in 0..8 {
        let mut applied = 0.;
        let mut applied_m = 0.;
        for load in loads {
            let (q, s) = match *load {
                FloorLoad::Area { slot: k, w, region } if k == slot => {
                    let [s0, s1, _, _] = region.unwrap_or([b.s0, b.s1, b.t0, b.t1]);
                    let q = w * (s1 - s0) * (region.map_or(b.t1 - b.t0, |x| x[3] - x[2]));
                    (q, (s0 + s1) * 0.5)
                }
                FloorLoad::Line {
                    slot: k,
                    w,
                    s,
                    t0,
                    t1,
                } if k == slot => (w * (t1 - t0), s),
                FloorLoad::LineAlong {
                    slot: k, w, s0, s1, ..
                } if k == slot => (w * (s1 - s0), (s0 + s1) * 0.5),
                FloorLoad::Point { slot: k, p, s, .. } if k == slot => (p, s),
                _ => continue,
            };
            applied += q;
            applied_m += q * s;
        }
        let mut resisted = 0.;
        let mut resisted_m = 0.;
        for (i, reaction) in r.per_support.iter().enumerate() {
            let q: f64 = reaction
                .lines
                .iter()
                .filter(|x| x.slot == slot)
                .map(|x| x.w * (x.t1 - x.t0))
                .sum::<f64>()
                + reaction
                    .points
                    .iter()
                    .filter(|x| x.slot == slot)
                    .map(|x| x.p)
                    .sum::<f64>();
            resisted += q;
            resisted_m += q * supports[i].s;
        }
        assert!(
            (resisted - applied).abs() <= 1e-12 * applied.abs().max(1.),
            "slot {slot}: {resisted} != {applied}"
        );
        assert!(
            (resisted_m - applied_m).abs() <= 1e-12 * applied_m.abs().max(1.),
            "slot {slot}: M {resisted_m} != {applied_m}"
        );
    }
}

/// Two simple spans under uniform pressure give half-span end reactions and the sum of adjacent half-spans, by hand.
#[test]
fn simple_two_span_uniform() {
    let r = transfer_floor(
        bounds(),
        &supports(),
        &[FloorLoad::Area {
            slot: 0,
            w: 2.,
            region: None,
        }],
        FloorContinuity::Simple,
    )
    .unwrap();
    equilibrium(
        bounds(),
        &supports(),
        &[FloorLoad::Area {
            slot: 0,
            w: 2.,
            region: None,
        }],
        &r,
    );
    for (i, want) in [16., 40., 24.].into_iter().enumerate() {
        let got: f64 = r.per_support[i]
            .lines
            .iter()
            .map(|x| x.w * (x.t1 - x.t0))
            .sum();
        assert!((got - want).abs() < 1e-12);
    }
    assert!((sum(&r, 0) - 80.).abs() < 1e-12);
}

/// Two equal continuous spans give the classical 3/8, 10/8, 3/8 reactions under uniform load.
#[test]
fn continuous_two_span_uniform() {
    let s = [
        FloorSupport {
            s: 0.,
            t0: 0.,
            t1: 4.,
        },
        FloorSupport {
            s: 5.,
            t0: 0.,
            t1: 4.,
        },
        FloorSupport {
            s: 10.,
            t0: 0.,
            t1: 4.,
        },
    ];
    let r = transfer_floor(
        bounds(),
        &s,
        &[FloorLoad::Area {
            slot: 0,
            w: 1.,
            region: Some([0., 10., 0., 1.]),
        }],
        FloorContinuity::Continuous,
    )
    .unwrap();
    equilibrium(
        FloorBounds {
            s0: 0.,
            s1: 10.,
            t0: 0.,
            t1: 4.,
        },
        &s,
        &[FloorLoad::Area {
            slot: 0,
            w: 1.,
            region: Some([0., 10., 0., 1.]),
        }],
        &r,
    );
    let got: Vec<f64> = r
        .per_support
        .iter()
        .map(|s| s.lines.iter().map(|x| x.w * (x.t1 - x.t0)).sum::<f64>())
        .collect();
    for (g, w) in got.into_iter().zip([1.875, 6.25, 1.875]) {
        assert!((g - w).abs() < 1e-11, "{g} != {w}");
    }
}

/// Unequal continuous spans match the three-moment equation for uniform load.
#[test]
fn continuous_unequal_spans_three_moment() {
    let s = [
        FloorSupport {
            s: 0.,
            t0: 0.,
            t1: 1.,
        },
        FloorSupport {
            s: 4.,
            t0: 0.,
            t1: 1.,
        },
        FloorSupport {
            s: 10.,
            t0: 0.,
            t1: 1.,
        },
    ];
    let r = transfer_floor(
        FloorBounds {
            s0: 0.,
            s1: 10.,
            t0: 0.,
            t1: 1.,
        },
        &s,
        &[FloorLoad::Area {
            slot: 0,
            w: 1.,
            region: None,
        }],
        FloorContinuity::Continuous,
    )
    .unwrap();
    equilibrium(
        FloorBounds {
            s0: 0.,
            s1: 10.,
            t0: 0.,
            t1: 1.,
        },
        &s,
        &[FloorLoad::Area {
            slot: 0,
            w: 1.,
            region: None,
        }],
        &r,
    );
    let m = -(4f64.powi(3) + 6f64.powi(3)) / (8. * (4. + 6.));
    let expected = [
        2. + m / 4.,
        10. - (2. + m / 4.) - (3. + m / 6.),
        3. + m / 6.,
    ];
    for (reaction, want) in r.per_support.iter().zip(expected) {
        let got: f64 = reaction.lines.iter().map(|x| x.w * (x.t1 - x.t0)).sum();
        assert!((got - want).abs() < 1e-11, "{got} != {want}");
    }
}

/// A balcony overhang relieves the back support; for two supports the continuous solution is the same determinate beam.
#[test]
fn balcony_overhang_has_back_uplift() {
    let b = FloorBounds {
        s0: 0.,
        s1: 13.,
        t0: 0.,
        t1: 1.,
    };
    let s = [
        FloorSupport {
            s: 0.,
            t0: 0.,
            t1: 1.,
        },
        FloorSupport {
            s: 5.,
            t0: 0.,
            t1: 1.,
        },
    ];
    for continuity in [FloorContinuity::Simple, FloorContinuity::Continuous] {
        let r = transfer_floor(
            b,
            &s,
            &[FloorLoad::Area {
                slot: 0,
                w: 1.,
                region: None,
            }],
            continuity,
        )
        .unwrap();
        equilibrium(
            b,
            &s,
            &[FloorLoad::Area {
                slot: 0,
                w: 1.,
                region: None,
            }],
            &r,
        );
        let reactions: Vec<f64> = r
            .per_support
            .iter()
            .map(|x| x.lines.iter().map(|l| l.w).sum())
            .collect();
        assert!((reactions[0] + 3.9).abs() < 1e-11);
        assert!((reactions[1] - 16.9).abs() < 1e-11);
        assert!((reactions.iter().sum::<f64>() - 13.).abs() < 1e-11);
    }
}

/// Point and joist-line loads use lever-rule reactions at their exact t coordinate.
#[test]
fn concentrated_loads_use_exact_t() {
    let r = transfer_floor(
        bounds(),
        &supports()[..2],
        &[
            FloorLoad::Point {
                slot: 0,
                p: 10.,
                s: 2.,
                t: 1.,
            },
            FloorLoad::LineAlong {
                slot: 1,
                w: 3.,
                t: 3.,
                s0: 1.,
                s1: 3.,
            },
            FloorLoad::Line {
                slot: 2,
                w: 2.,
                s: 3.,
                t0: 0.,
                t1: 4.,
            },
        ],
        FloorContinuity::Simple,
    )
    .unwrap();
    equilibrium(
        bounds(),
        &supports()[..2],
        &[
            FloorLoad::Point {
                slot: 0,
                p: 10.,
                s: 2.,
                t: 1.,
            },
            FloorLoad::LineAlong {
                slot: 1,
                w: 3.,
                t: 3.,
                s0: 1.,
                s1: 3.,
            },
            FloorLoad::Line {
                slot: 2,
                w: 2.,
                s: 3.,
                t0: 0.,
                t1: 4.,
            },
        ],
        &r,
    );
    assert_eq!(
        r.per_support[0]
            .points
            .iter()
            .find(|p| p.slot == 0)
            .unwrap()
            .t,
        1.
    );
    let p0: f64 = r.per_support[0]
        .points
        .iter()
        .filter(|p| p.slot == 0)
        .map(|p| p.p)
        .sum();
    let p1: f64 = r.per_support[1]
        .points
        .iter()
        .filter(|p| p.slot == 0)
        .map(|p| p.p)
        .sum();
    assert!((p0 - 5.).abs() < 1e-12 && (p1 - 5.).abs() < 1e-12);
    assert!((sum(&r, 1) - 6.).abs() < 1e-12);
    assert!((sum(&r, 2) - 8.).abs() < 1e-12);
}

/// A partial bearing is split into independent strips, and each reaction maps to wall-local downward load.
#[test]
fn partial_support_and_wall_mapping() {
    let s = [
        FloorSupport {
            s: 0.,
            t0: 0.,
            t1: 4.,
        },
        FloorSupport {
            s: 5.,
            t0: 0.,
            t1: 2.,
        },
        FloorSupport {
            s: 10.,
            t0: 0.,
            t1: 4.,
        },
    ];
    let r = transfer_floor(
        bounds(),
        &s,
        &[FloorLoad::Area {
            slot: 0,
            w: 1.,
            region: None,
        }],
        FloorContinuity::Simple,
    )
    .unwrap();
    equilibrium(
        bounds(),
        &s,
        &[FloorLoad::Area {
            slot: 0,
            w: 1.,
            region: None,
        }],
        &r,
    );
    assert!((sum(&r, 0) - 40.).abs() < 1e-12);
    let loads = r.per_support[0].wall_loads(100., true);
    assert!(loads[0]
        .iter()
        .any(|l| matches!(l,WallLoad::TopLine{wy0,..} if *wy0<0.)));
}

/// Unsupported loaded strips, a single support, out-of-floor loads, and invalid slots are errors.
#[test]
fn invalid_floor_inputs_error() {
    let unsupported = [FloorSupport {
        s: 0.,
        t0: 0.,
        t1: 2.,
    }];
    assert!(transfer_floor(
        bounds(),
        &unsupported,
        &[FloorLoad::Area {
            slot: 0,
            w: 1.,
            region: None
        }],
        FloorContinuity::Simple
    )
    .is_err());
    assert!(transfer_floor(
        bounds(),
        &supports()[..1],
        &[FloorLoad::Point {
            slot: 0,
            p: 1.,
            s: 1.,
            t: 1.
        }],
        FloorContinuity::Simple
    )
    .is_err());
    assert!(transfer_floor(
        bounds(),
        &supports(),
        &[FloorLoad::Point {
            slot: 0,
            p: 1.,
            s: 11.,
            t: 1.
        }],
        FloorContinuity::Simple
    )
    .is_err());
    assert!(transfer_floor(
        bounds(),
        &supports(),
        &[FloorLoad::Point {
            slot: 8,
            p: 1.,
            s: 1.,
            t: 1.
        }],
        FloorContinuity::Simple
    )
    .is_err());
}

/// Floor slot reactions mapped to a window wall give the hand header takedown, including a factored D/L combination (1e-9).
#[test]
fn floor_to_window_wall_factored_takedown() {
    let (length, height) = (3600., 2440.);
    let support = [
        FloorSupport {
            s: 0.,
            t0: 0.,
            t1: length,
        },
        FloorSupport {
            s: 3000.,
            t0: 0.,
            t1: length,
        },
    ];
    let reactions = transfer_floor(
        FloorBounds {
            s0: 0.,
            s1: 3000.,
            t0: 0.,
            t1: length,
        },
        &support,
        &[
            FloorLoad::Area {
                slot: 0,
                w: 10. / 1500.,
                region: None,
            },
            FloorLoad::Area {
                slot: 1,
                w: 4. / 1500.,
                region: None,
            },
        ],
        FloorContinuity::Simple,
    )
    .unwrap();
    equilibrium(
        FloorBounds {
            s0: 0.,
            s1: 3000.,
            t0: 0.,
            t1: length,
        },
        &support,
        &[
            FloorLoad::Area {
                slot: 0,
                w: 10. / 1500.,
                region: None,
            },
            FloorLoad::Area {
                slot: 1,
                w: 4. / 1500.,
                region: None,
            },
        ],
        &reactions,
    );
    let mut wall_loads = reactions.per_support[0].wall_loads(0., false);
    let studs = [
        0., 400., 800., 1100., 1200., 1600., 2000., 2400., 2500., 2800., 3200., 3600.,
    ]
    .into_iter()
    .map(|x| WallStud { x, ea: 50e6 })
    .collect();
    let layout = StudWallLayout {
        length,
        height,
        studs,
        top_ea: 100e6,
        top_ei: None,
        top_distribution: ChordDistribution::Direct,
        bottom_ea: 50e6,
        bottom_ei: 6e9,
        openings: vec![WallOpening {
            x0: 1140.,
            x1: 2460.,
            head_height: 2000.,
            header_ei: 7.8e11,
            jack_ea: 50e6,
            door: false,
        }],
        sheathing: None,
        foundation: true,
        floor: None,
    };
    let mut d: Domain<2, 3, NodeId, Element> = Domain::new();
    let corners = [
        d.add_node(Node::new([0., 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, 0.]).fix(0).fix(1)),
        d.add_node(Node::new([length, height])),
        d.add_node(Node::new([0., height]).fix(0)),
    ];
    let wall = d.add_element(Element::StudWall(
        StudWall::new(corners, vec![], layout, std::mem::take(&mut wall_loads)).unwrap(),
    ));
    let pattern = d.default_pattern();
    let mut factors = [0.; WALL_SLOTS];
    factors[0] = 1.2;
    factors[1] = 1.6;
    d.add_element_load(
        pattern,
        wall,
        ElementLoad {
            wall: factors,
            ..Default::default()
        },
    );
    let mut a = AnalysisBuilder::new()
        .constraint_handler(ConstraintHandler::Transformation)
        .integrator(Integrator::LoadControl { increment: 1. })
        .algorithm(Algorithm::Linear)
        .test(ConvergenceTest::NormUnbalance {
            tol: 1e-11,
            max_iter: 4,
        })
        .build(d);
    a.step().unwrap();
    let response = a.domain().element_stud_wall_response(wall, 1.).unwrap();
    for jack in &response.jacks {
        assert!(
            (jack.axial + 11_960.).abs() <= 11_960. * 1e-9,
            "jack at {} = {}",
            jack.x,
            jack.axial
        );
    }
    let base: f64 = response.base_nodes.iter().map(|r| r.fy).sum();
    assert!((base - 18.4 * length).abs() <= 18.4 * length * 1e-9);
}
