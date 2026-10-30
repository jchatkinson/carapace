//! Section probe: analytic cases on small sections. JS-engine parity lives in `section_probe_parity.rs`.

use carapace_wasm::input_v1::materials::MaterialSpec;
use carapace_wasm::section_probe::*;

/// `n` fibers across a rectangle of `width` x `depth` (y is the depth direction), single material `0`.
fn rectangle(width: f64, depth: f64, n: usize) -> ProbeFibers {
    let dy = depth / n as f64;
    ProbeFibers {
        y: (0..n)
            .map(|i| -depth / 2.0 + dy * (i as f64 + 0.5))
            .collect(),
        z: vec![0.0; n],
        area: vec![width * dy; n],
        material: vec![0; n],
        group: vec![],
    }
}

fn elastic_probe(e: f64, n: usize) -> SectionProbe {
    SectionProbe::new(&SectionProbeConfig {
        core_materials: vec![MaterialSpec::Elastic { e }],
        materials: vec![ProbeMaterialSpec::Core { index: 0 }],
        fibers: rectangle(10.0, 20.0, n),
    })
    .unwrap()
}

fn near(got: f64, want: f64, rel: f64) {
    let scale = want.abs().max(1.0);
    assert!((got - want).abs() <= rel * scale, "{got} vs {want}");
}

/// Elastic rectangle: axial force is EA*eps0 and the bending moment about z is EI*kz (sign: mz = -sum(y f)).
#[test]
fn elastic_resultants_match_closed_form() {
    let p = elastic_probe(30_000.0, 400);
    let (a, i) = (10.0 * 20.0, 10.0 * 20.0_f64.powi(3) / 12.0);
    let r = p.evaluate(1e-4, 0.0, 0.0);
    near(r.n, 30_000.0 * a * 1e-4, 1e-9);
    near(r.mz, 0.0, 1e-9);
    let r = p.evaluate(0.0, 2e-5, 0.0);
    near(r.n, 0.0, 1e-9);
    // strain = -y*kz, force = -E y kz dA, mz = -sum(y f) = +E I kz; midpoint rule is exact for this linear case
    near(r.mz, 30_000.0 * i * 2e-5, 1e-3);
    near(r.m2, -r.mz, 1e-12);
    assert_eq!(r.m3, r.my);
}

/// Per-group sums add up to the total, and the group stride is `3 * (1 + groups)`.
#[test]
fn group_sums_add_to_total() {
    let mut fibers = rectangle(10.0, 20.0, 40);
    fibers.group = (0..40).map(|i| (i % 3) as u32).collect();
    let p = SectionProbe::new(&SectionProbeConfig {
        core_materials: vec![MaterialSpec::Elastic { e: 1000.0 }],
        materials: vec![ProbeMaterialSpec::Core { index: 0 }],
        fibers,
    })
    .unwrap();
    assert_eq!(p.group_count(), 3);
    let out = p.evaluate_batch(&[1e-3, 1e-5, 2e-5], true);
    assert_eq!(out.len(), 3 * (1 + 3));
    for c in 0..3 {
        near(out[c], out[3 + c] + out[6 + c] + out[9 + c], 1e-12);
    }
    let r = p.evaluate(1e-3, 1e-5, 2e-5);
    near(out[0], r.n, 1e-12);
    near(out[1], r.mz, 1e-12);
    near(out[2], r.my, 1e-12);
}

/// The probe-only curves reproduce hand-computed values, including the Whitney step, strain cutoffs, user
/// curves with and without a drop outside the table, and the post-crush descent.
#[test]
fn probe_only_materials() {
    let build = |m: ProbeMaterialSpec| {
        SectionProbe::new(&SectionProbeConfig {
            core_materials: vec![],
            materials: vec![m],
            fibers: ProbeFibers {
                y: vec![0.0],
                z: vec![0.0],
                area: vec![1.0],
                material: vec![0],
                group: vec![],
            },
        })
        .unwrap()
    };
    let stress = |p: &SectionProbe, e: f64| p.evaluate(e, 0.0, 0.0).n;
    let concrete = |model| ConcreteSpec {
        model,
        fc: 4000.0,
        ec: 3_600_000.0,
        ecu: 0.003,
        eu: Some(0.006),
        fu: 400.0,
        alpha: 0.85,
        beta1: Some(0.85),
        tension_capacity: false,
        ft: 0.0,
        ets: 0.0,
        curve_points: vec![],
        drop_after_last_point: false,
    };
    // Whitney: zero below ecu*(1-beta1) = 0.00045, alpha*fc beyond, linear descent to -fu between ecu and eu
    let w = build(ProbeMaterialSpec::Concrete(concrete(
        ConcreteModel::Whitney,
    )));
    near(stress(&w, -0.0004), 0.0, 1e-12);
    near(stress(&w, -0.0005), -0.85 * 4000.0, 1e-12);
    near(stress(&w, 0.001), 0.0, 1e-12);
    // halfway between ecu and eu the stress is halfway between -alpha fc and -fu
    near(stress(&w, -0.0045), 0.5 * (-0.85 * 4000.0 - 400.0), 1e-9);
    near(stress(&w, -0.01), -400.0, 1e-12);
    // Hognestad: peak fc at e0 = 2 fc / Ec, 0.85 fc at ecu
    let e0 = 2.0 * 4000.0 / 3_600_000.0;
    let h = build(ProbeMaterialSpec::Concrete(concrete(
        ConcreteModel::Hognestad,
    )));
    near(stress(&h, -e0), -4000.0, 1e-12);
    near(stress(&h, -0.003), -0.85 * 4000.0, 1e-9);
    // Parabolic holds fc beyond e0
    let pb = build(ProbeMaterialSpec::Concrete(concrete(
        ConcreteModel::Parabolic,
    )));
    near(stress(&pb, -0.0025), -4000.0, 1e-12);
    // tension branch: linear to ft, then softening at ets
    let mut t = concrete(ConcreteModel::Hognestad);
    (t.tension_capacity, t.ft, t.ets) = (true, 400.0, 100_000.0);
    let t = build(ProbeMaterialSpec::Concrete(t));
    near(stress(&t, 0.0001), 360.0, 1e-9);
    near(
        stress(&t, 400.0 / 3_600_000.0 + 0.002),
        400.0 - 100_000.0 * 0.002,
        1e-6,
    );
    near(stress(&t, 0.1), 0.0, 1e-12);
    // steel: hardening 2% and a tensile cutoff
    let steel = |hardening, eps_ut| {
        build(ProbeMaterialSpec::Steel(SteelSpec {
            model: SteelModel::ElasticPlastic,
            fy: 60_000.0,
            es: 30_000_000.0,
            hardening_ratio: hardening,
            eps_uc: 0.0,
            eps_ut,
            curve_points: vec![],
            drop_after_last_point: false,
        }))
    };
    let s = steel(0.02, 0.05);
    near(stress(&s, 0.001), 30_000.0, 1e-12);
    near(
        stress(&s, 0.004),
        60_000.0 + 0.02 * 30_000_000.0 * 0.002,
        1e-12,
    );
    near(
        stress(&s, -0.004),
        -(60_000.0 + 0.02 * 30_000_000.0 * 0.002),
        1e-12,
    );
    near(stress(&s, 0.06), 0.0, 1e-12);
    // user curve: interpolates, holds the end value, or drops to zero outside the table
    let curve = |drop| {
        build(ProbeMaterialSpec::General(GeneralSpec {
            model: GeneralModel::UserCurve,
            e: 0.0,
            eps_uc: 0.0,
            eps_ut: 0.0,
            curve_points: vec![[0.0, 0.0], [0.01, 100.0], [0.005, 40.0]],
            drop_after_last_point: drop,
        }))
    };
    let held = curve(false);
    near(stress(&held, 0.0075), 70.0, 1e-12); // points are sorted: (0,0), (0.005,40), (0.01,100)
    near(stress(&held, 0.02), 100.0, 1e-12);
    near(stress(&held, -1.0), 0.0, 1e-12);
    let dropped = curve(true);
    near(stress(&dropped, 0.02), 0.0, 1e-12);
    near(stress(&dropped, -0.01), 0.0, 1e-12);
}

/// A Whitney concrete without beta1 and a combine term that is not earlier are rejected, as is a bad core index.
#[test]
fn invalid_materials_are_rejected() {
    let fibers = || ProbeFibers {
        y: vec![0.0],
        z: vec![0.0],
        area: vec![1.0],
        material: vec![0],
        group: vec![],
    };
    let cfg = |materials| SectionProbeConfig {
        core_materials: vec![],
        materials,
        fibers: fibers(),
    };
    let whitney = ConcreteSpec {
        model: ConcreteModel::Whitney,
        fc: 4000.0,
        ec: 3.6e6,
        ecu: 0.003,
        eu: None,
        fu: 0.0,
        alpha: 0.85,
        beta1: None,
        tension_capacity: false,
        ft: 0.0,
        ets: 0.0,
        curve_points: vec![],
        drop_after_last_point: false,
    };
    assert!(SectionProbe::new(&cfg(vec![ProbeMaterialSpec::Concrete(whitney)])).is_err());
    assert!(SectionProbe::new(&cfg(vec![ProbeMaterialSpec::Combine {
        terms: vec![(0, 1.0)]
    }]))
    .is_err());
    assert!(SectionProbe::new(&cfg(vec![ProbeMaterialSpec::Core { index: 0 }])).is_err());
}

/// A bar fiber defined as steel minus the displaced concrete carries (Es - Ec) * eps * A while both are elastic.
#[test]
fn combine_subtracts_the_host() {
    let p = SectionProbe::new(&SectionProbeConfig {
        core_materials: vec![
            MaterialSpec::Elastic { e: 30_000.0 },
            MaterialSpec::Elastic { e: 29_000_000.0 },
        ],
        materials: vec![
            ProbeMaterialSpec::Core { index: 0 },
            ProbeMaterialSpec::Core { index: 1 },
            ProbeMaterialSpec::Combine {
                terms: vec![(1, 1.0), (0, -1.0)],
            },
        ],
        fibers: ProbeFibers {
            y: vec![0.0],
            z: vec![0.0],
            area: vec![2.0],
            material: vec![2],
            group: vec![],
        },
    })
    .unwrap();
    near(
        p.evaluate(0.001, 0.0, 0.0).n,
        (29_000_000.0 - 30_000.0) * 0.001 * 2.0,
        1e-12,
    );
}

/// A profile (et, eb) at angle theta puts et at the fibre with the largest r and eb at the smallest.
#[test]
fn profile_states_hit_the_extreme_strains() {
    let p = elastic_probe(1.0, 50);
    let y: Vec<f64> = (0..50).map(|i| -10.0 + 0.4 * (i as f64 + 0.5)).collect();
    for theta in [0.0, 0.7, 1.9, std::f64::consts::PI] {
        let (e0, kz, ky) = p.profile_state(theta, -0.0035, 0.002).unwrap();
        let strains = p.fiber_responses(e0, kz, ky);
        let r: Vec<f64> = y.iter().map(|y| y * theta.cos()).collect();
        if theta.cos().abs() < 1e-6 {
            continue; // all fibres on the y axis share r = 0
        }
        let (imax, imin) = (
            (0..50)
                .max_by(|a, b| r[*a].partial_cmp(&r[*b]).unwrap())
                .unwrap(),
            (0..50)
                .min_by(|a, b| r[*a].partial_cmp(&r[*b]).unwrap())
                .unwrap(),
        );
        near(strains[2 * imax], -0.0035, 1e-9);
        near(strains[2 * imin], 0.002, 1e-9);
    }
}

/// Uniform-strain axial extremes of an elastic-perfectly-plastic steel section are +-fy*A at yield.
#[test]
fn axial_extremes_of_an_elastic_plastic_section() {
    let p = SectionProbe::new(&SectionProbeConfig {
        core_materials: vec![],
        materials: vec![ProbeMaterialSpec::Steel(SteelSpec {
            model: SteelModel::ElasticPlastic,
            fy: 50.0,
            es: 29_000.0,
            hardening_ratio: 0.0,
            eps_uc: 0.0,
            eps_ut: 0.0,
            curve_points: vec![],
            drop_after_last_point: false,
        })],
        fibers: rectangle(10.0, 20.0, 20),
    })
    .unwrap();
    let ext = p.axial_extremes(0.02).unwrap();
    near(ext.compression.response.n, -50.0 * 200.0, 1e-9);
    near(ext.tension.response.n, 50.0 * 200.0, 1e-9);
}

/// Elastic moment-curvature at constant axial force: eps0 = P / EA and M = EI * kappa.
#[test]
fn elastic_moment_curvature() {
    let p = elastic_probe(30_000.0, 200);
    let (a, i) = (200.0, 10.0 * 8000.0 / 12.0);
    let res = p
        .moment_curvature(&MomentCurvatureRequest {
            axial: -5_000.0,
            angle: 0.0,
            curvatures: vec![0.0, 1e-5, 2e-5, 4e-5],
            tolerance: 1e-9,
        })
        .unwrap();
    assert_eq!(res.non_converged, 0);
    for pt in &res.points {
        near(pt.eps0, -5_000.0 / (30_000.0 * a), 1e-6);
        near(pt.m2.abs(), 30_000.0 * i * pt.curvature, 1e-3);
        assert!(pt.min_strain <= pt.max_strain);
    }
}

/// Every solved isoline point satisfies N = P and has its moment along the requested angle; points
/// with the section's axial range inside a level are found, and the outer root is returned.
#[test]
fn isolines_satisfy_their_equations() {
    // symmetric rectangle, elastic-plastic steel-like material, tension and compression capable
    let p = SectionProbe::new(&SectionProbeConfig {
        core_materials: vec![],
        materials: vec![ProbeMaterialSpec::Steel(SteelSpec {
            model: SteelModel::ElasticPlastic,
            fy: 50.0,
            es: 29_000.0,
            hardening_ratio: 0.0,
            eps_uc: 0.0,
            eps_ut: 0.0,
            curve_points: vec![],
            drop_after_last_point: false,
        })],
        fibers: ProbeFibers {
            // 6 x 10 grid of fibers
            y: (0..60)
                .map(|k| -5.0 + 10.0 * ((k / 6) as f64 + 0.5) / 10.0)
                .collect(),
            z: (0..60)
                .map(|k| -3.0 + 6.0 * ((k % 6) as f64 + 0.5) / 6.0)
                .collect(),
            area: vec![1.0; 60],
            material: vec![0; 60],
            group: vec![],
        },
    })
    .unwrap();
    // profile loop for ecu = 0.003, steel strain limit 0.03: compression face at ecu, far face from ecu to +0.03,
    // then the mirrored half, so the loop is symmetric under (et, eb) -> (eb, et) reversed.
    let ecu = -0.003;
    let far = [ecu, -0.002, -0.001, 0.0, 0.002, 0.005, 0.01, 0.02, 0.03];
    let mut profiles: Vec<f64> = Vec::new();
    for e in far {
        profiles.extend_from_slice(&[ecu, e]);
    }
    for e in far.iter().rev() {
        profiles.extend_from_slice(&[*e, ecu]);
    }
    let angles: Vec<f64> = (0..36)
        .map(|i| std::f64::consts::PI * i as f64 / 36.0)
        .collect();
    let mut p = p;
    let count = p.build_grid(&angles, &profiles, None).unwrap();
    assert_eq!(count, 36 * 18);
    let (lo, hi) = p.grid_axial_range(false).unwrap();
    let levels: Vec<f64> = (1..9).map(|k| lo + (hi - lo) * k as f64 / 9.0).collect();
    let moment_angles = [0.0, 0.4, 1.1, 2.0];
    let res = p
        .solve_isolines(&IsolineRequest {
            levels: levels.clone(),
            moment_angles: moment_angles.to_vec(),
            options: IsolineOptions::default(),
            factored: false,
        })
        .unwrap();
    assert_eq!(res.len(), levels.len() * moment_angles.len());
    let mut found = 0;
    for pt in &res {
        if !pt.found {
            continue;
        }
        found += 1;
        let target = levels[pt.level_index as usize];
        let alpha = moment_angles[pt.moment_angle_index as usize];
        near(pt.n, target, 1e-4 * (hi - lo));
        // moment direction: cross product with (cos a, sin a) is ~0 and the dot product is positive
        let cross = pt.m2 * alpha.sin() - pt.m3 * alpha.cos();
        assert!(cross.abs() < 1e-3 * pt.moment.max(1.0), "{pt:?}");
        assert!(pt.moment >= 0.0, "{pt:?}");
        // the state really is at the crushing limit: one of the extreme fibres is at ecu
        assert!(
            (pt.top_strain - ecu).abs() < 1e-9 || (pt.bottom_strain - ecu).abs() < 1e-9,
            "{pt:?}"
        );
    }
    assert!(
        found >= res.len() - 2,
        "only {found} of {} solved",
        res.len()
    );
}

/// 6 x 10 elastic-plastic grid with the outer rows in group 1 ("bars") and the rest in group 0.
fn grouped_grid() -> SectionProbe {
    let y: Vec<f64> = (0..60)
        .map(|k| -5.0 + 10.0 * ((k / 6) as f64 + 0.5) / 10.0)
        .collect();
    let group = (0..60)
        .map(|k| u32::from(k / 6 == 0 || k / 6 == 9))
        .collect();
    SectionProbe::new(&SectionProbeConfig {
        core_materials: vec![],
        materials: vec![ProbeMaterialSpec::Steel(SteelSpec {
            model: SteelModel::ElasticPlastic,
            fy: 50.0,
            es: 29_000.0,
            hardening_ratio: 0.0,
            eps_uc: 0.0,
            eps_ut: 0.0,
            curve_points: vec![],
            drop_after_last_point: false,
        })],
        fibers: ProbeFibers {
            y,
            z: (0..60)
                .map(|k| -3.0 + 6.0 * ((k % 6) as f64 + 0.5) / 6.0)
                .collect(),
            area: vec![1.0; 60],
            material: vec![0; 60],
            group,
        },
    })
    .unwrap()
}

/// Group factors scale each group's force; ACI phi follows the bar group's largest tensile strain.
#[test]
fn factoring_scales_groups_and_follows_the_aci_ramp() {
    let mut p = grouped_grid();
    let (e0, kz, ky) = (0.0004, 0.0, 0.0);
    let mut g = [0.0; 9];
    p.evaluate_grouped(e0, kz, ky, &mut g);
    let by_group = p.evaluate_factored(
        &Factoring::Groups {
            factors: vec![0.65, 0.85],
        },
        e0,
        kz,
        ky,
    );
    near(by_group.response.n, 0.65 * g[3] + 0.85 * g[6], 1e-12);
    near(by_group.response.mz, 0.65 * g[4] + 0.85 * g[7], 1e-12);
    assert!(by_group.phi.is_none());

    let aci = |eps: f64| {
        p.evaluate_factored(
            &Factoring::Aci {
                bar_group: 1,
                phi_tension: 0.9,
                phi_compression: 0.65,
                eps_y: 0.002,
                eps_tl: 0.005,
            },
            eps,
            0.0,
            0.0,
        )
    };
    // uniform strain: every bar has strain eps
    assert_eq!(aci(0.0015).phi, Some(0.65));
    assert_eq!(aci(-0.002).phi, Some(0.65));
    assert_eq!(aci(-0.002).eps_t, Some(0.0));
    near(aci(0.0035).phi.unwrap(), 0.775, 1e-12); // halfway between eps_y and eps_tl
    assert_eq!(aci(0.01).phi, Some(0.9));
    let nominal = p.evaluate(0.0035, 0.0, 0.0);
    near(aci(0.0035).response.n, 0.775 * nominal.n, 1e-12);
    assert!(p
        .build_grid(
            &[0.0],
            &[0.0; 6],
            Some(Factoring::Aci {
                bar_group: 5,
                phi_tension: 0.9,
                phi_compression: 0.65,
                eps_y: 0.002,
                eps_tl: 0.005
            })
        )
        .is_err());
}

/// Factored isolines hit the requested factored axial force, and sample_states reports phi and c.
#[test]
fn factored_isolines_solve_in_factored_space() {
    let mut p = grouped_grid();
    let ecu = -0.003;
    let far = [ecu, -0.002, -0.001, 0.0, 0.002, 0.005, 0.01, 0.02, 0.03];
    let mut profiles: Vec<f64> = Vec::new();
    for e in far {
        profiles.extend_from_slice(&[ecu, e]);
    }
    for e in far.iter().rev() {
        profiles.extend_from_slice(&[*e, ecu]);
    }
    let angles: Vec<f64> = (0..36)
        .map(|i| std::f64::consts::PI * i as f64 / 36.0)
        .collect();
    let aci = Factoring::Aci {
        bar_group: 1,
        phi_tension: 0.9,
        phi_compression: 0.65,
        eps_y: 50.0 / 29_000.0,
        eps_tl: 0.005,
    };
    let rows = p.sample_states(&angles, &profiles, Some(&aci)).unwrap();
    assert_eq!(rows.len(), 36 * 18 * STATE_STRIDE);
    let row = &rows[STATE_STRIDE * 5..STATE_STRIDE * 6];
    let phi = row[6];
    assert!((0.65..=0.9).contains(&phi), "{phi}");
    near(row[3], phi * row[0], 1e-12);
    near(row[4], phi * row[1], 1e-12);
    // first profile is uniform crushing: c undefined; the second has the neutral axis outside the section
    assert!(rows[8].is_nan());
    assert_eq!((rows[9], rows[10]), (ecu, ecu));

    p.build_grid(&angles, &profiles, Some(aci)).unwrap();
    let (lo, hi) = p.grid_axial_range(true).unwrap();
    let (nlo, nhi) = p.grid_axial_range(false).unwrap();
    assert!(
        hi <= nhi + 1e-9 && lo >= nlo - 1e-9,
        "{lo} {hi} {nlo} {nhi}"
    );
    let levels: Vec<f64> = (1..9).map(|k| lo + (hi - lo) * k as f64 / 9.0).collect();
    let res = p
        .solve_isolines(&IsolineRequest {
            levels: levels.clone(),
            moment_angles: vec![0.0, 1.0],
            options: IsolineOptions::default(),
            factored: true,
        })
        .unwrap();
    let mut found = 0;
    for pt in res.iter().filter(|pt| pt.found) {
        found += 1;
        near(pt.n, levels[pt.level_index as usize], 1e-4 * (hi - lo));
        let phi = pt.phi.expect("ACI phi is reported for factored points");
        assert!((0.65 - 1e-12..=0.9 + 1e-12).contains(&phi), "{pt:?}");
        assert!(pt.eps_t.is_some());
    }
    assert!(found >= res.len() * 3 / 4, "{found} of {}", res.len());
    // factored isolines need a factoring
    let mut q = grouped_grid();
    q.build_grid(&angles, &profiles, None).unwrap();
    assert!(q
        .solve_isolines(&IsolineRequest {
            levels: vec![0.0],
            moment_angles: vec![0.0],
            options: IsolineOptions::default(),
            factored: true,
        })
        .is_err());
}
