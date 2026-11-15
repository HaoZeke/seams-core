/* Copyright (c) d-SEAMS developers */
fn diamond(nrep: usize) -> (Vec<f64>, f64) {
    let a = 2.75 * 4.0 / 3.0_f64.sqrt();
    let fcc = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]];
    let basis = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]];
    let mut xyz = Vec::new();
    for i in 0..nrep {
        for j in 0..nrep {
            for k in 0..nrep {
                for shift in fcc {
                    for site in basis {
                        xyz.push((i as f64 + shift[0] + site[0]) * a);
                        xyz.push((j as f64 + shift[1] + site[1]) * a);
                        xyz.push((k as f64 + shift[2] + site[2]) * a);
                    }
                }
            }
        }
    }
    (xyz, nrep as f64 * a)
}

#[test]
fn cubic_ice_through_the_c_api() {
    seams::set_num_threads(1).unwrap();
    assert_eq!(seams::num_threads(), 1);
    assert!(!seams::version().is_empty());
    let (xyz, length) = diamond(2);
    let mut sys = seams::System::new().unwrap();
    sys.set_ortho(&xyz, length, length, length).unwrap();
    let (labels, counts) = sys
        .chill_plus(3.5, seams_sys::SEAMS_GRAPH_CUTOFF, seams_sys::SEAMS_CHILL_PUBLISHED)
        .unwrap();
    assert!(labels.iter().all(|label| *label == 0));
    assert_eq!(counts[0] as usize, labels.len());
    assert_eq!(sys.count() as usize, labels.len());
    let (_, counts) = sys
        .chill_plus(3.5, seams_sys::SEAMS_GRAPH_MUTUAL, seams_sys::SEAMS_CHILL_OVITO)
        .unwrap();
    assert_eq!(counts[0] as usize, labels.len());
    let (q4, q6) = sys.q4q6(3.5).unwrap();
    let mean = q6.iter().sum::<f64>() / q6.len() as f64;
    assert!(mean > 0.3);
    assert!(q4.iter().all(|v| v.is_finite()));
    let degree = sys.neighbour_degree(3.5).unwrap();
    assert!(degree.iter().all(|d| *d == 4));
    let (six, fraction) = sys.six_ring_fraction(3.5).unwrap();
    assert!(six > 0);
    assert!(fraction > 0.5);
    let mut tilted = [0.0; 9];
    tilted[3] = length;
    tilted[4] = length;
    tilted[5] = length;
    tilted[6] = 0.3;
    sys.set_frame(&xyz, &tilted).unwrap();
    sys.chill_plus(3.5, seams_sys::SEAMS_GRAPH_CUTOFF, seams_sys::SEAMS_CHILL_PUBLISHED)
        .unwrap();
}

fn chordless_reference(path: &str) -> (usize, usize, Vec<i32>, Vec<i32>) {
    let text = std::fs::read_to_string(path).unwrap();
    let nums: Vec<i32> = text.split_whitespace().map(|s| s.parse().unwrap()).collect();
    let (atoms, n, size) = (nums[0] as usize, nums[1] as usize, nums[2] as usize);
    let rings = nums[3..3 + n * size].to_vec();
    let per_atom = nums[3 + n * size..3 + n * size + atoms].to_vec();
    (atoms, size, rings, per_atom)
}

#[test]
fn chordless_six_cycles_of_cubic_ice() {
    let (xyz, length) = diamond(2);
    let mut sys = seams::System::new().unwrap();
    sys.set_ortho(&xyz, length, length, length).unwrap();
    let found = sys
        .rings_chordless(3.5, 6, 0, seams_sys::SEAMS_GRAPH_CUTOFF)
        .unwrap();
    assert_eq!(found.count(), 128);
    assert!(found.per_atom.iter().all(|n| *n == 12));
    for cycle in found.cycles() {
        assert!(cycle[1..].iter().all(|atom| *atom > cycle[0]));
    }
    assert_eq!(found.cycle(0), &found.rings[..6]);
    let nb = sys.neighbours(3.5, 0, seams_sys::SEAMS_GRAPH_CUTOFF).unwrap();
    assert_eq!(nb.chordless(4).unwrap().count(), 0);
    assert!(nb.chordless(13).is_err());
}

#[test]
fn chordless_matches_the_c_reference() {
    let (Ok(path), Ok(traj)) = (
        std::env::var("SEAMS_CHORDLESS_REF"),
        std::env::var("SEAMS_CHORDLESS_TRAJ"),
    ) else {
        eprintln!("SEAMS_CHORDLESS_REF is not set");
        return;
    };
    let (atoms, size, rings, per_atom) = chordless_reference(&path);
    let mut sys = seams::System::new().unwrap();
    sys.read(&traj, 1, 2).unwrap();
    assert_eq!(sys.count() as usize, atoms);
    let found = sys
        .rings_chordless(3.5, size as i32, 0, seams_sys::SEAMS_GRAPH_CUTOFF)
        .unwrap();
    assert_eq!(found.size, size);
    assert_eq!(found.rings, rings);
    assert_eq!(found.per_atom, per_atom);
}
