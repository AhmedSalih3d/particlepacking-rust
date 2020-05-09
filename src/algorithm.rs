use point;
use rayon::prelude::*;

#[derive(Copy, Clone)]
struct UpdatePoints{
    utmp: point::Point,
    ptmp: point::Point
}

struct AC;

// A struct with two fields
impl AC{
    // const H:f32  = 0.04;
    const H1:f32 = 1.0 / 0.04;
    const H2:f32 = 4.0 * 0.04 * 0.04;
    const AD:f32 = 348.154;
    const FAC:f32 = 5.0 / 8.0;
    const BETA:f32 = 4.0;
    const ZETA:f32 = 0.067_966;
    const V:f32 = 0.000_865_93;
    const DT:f32 = 0.014_713;
    const EPS:f32 = 0.000_001;
}

pub fn packstep_s(data: &mut point::Particles,
                  ncut: usize) {
    let n = data.pvec.len();

    // Indices ordered by x position
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| data.pvec[i].0.partial_cmp(&data.pvec[j].0).unwrap());
    // Reverse look-up of order
    let mut index_pos = vec![0; n];
    for i in 0..n {
        index_pos[order[i]] = i;
    }

    let mut pi_vec = data.pvec[0..(n-ncut)].to_vec();
    let mut ui_vec = data.uvec[0..(n-ncut)].to_vec();

    pi_vec
    .par_iter_mut()
    .zip(&mut ui_vec)
    .enumerate()
    .for_each(|(i,(p_ptr,u_ptr))|
        {
            let up_points = packstep_single(&data.pvec,&data.uvec,i,&order,&index_pos);
            *p_ptr = up_points.ptmp;
            *u_ptr = up_points.utmp;
        }
    );

    // Only alters relevant indices
    data.pvec[..(n - ncut)].clone_from_slice(&pi_vec);
    data.uvec[..(n - ncut)].clone_from_slice(&ui_vec);
}

fn packstep_single(pvec: &[point::Point], uvec: &[point::Point], iter: usize,  order: &[usize], index_pos: &[usize]) -> UpdatePoints{
    let mut wgx = 0.0;
    let mut wgy = 0.0;
    let mut wgz = 0.0;
    let p_i = pvec[iter];

    let mut close_points = Vec::new();

    let mut check = |j| -> bool {
        let p_j = &pvec[j];
        let rij = p_i - *p_j;
        let rij2 = (rij * rij).sum();
        if rij2 <= AC::H2 {
            close_points.push(j);
        }
        // If x distance is too large,
        // no more points are needed to checked
        rij.0 * rij.0 > AC::H2
    };

    let iter_pos = index_pos[iter];
    assert!(order[iter_pos] == iter);


    for &j in order[iter_pos+1..].iter() {
        if check(j) {
            break;
        }
    }

    for &j in order[..iter_pos].iter().rev() {
        if check(j) {
            break;
        }
    }

    // It is for if you want the exact floating point sum order (for validation)
    // close_points.sort();

    for j in close_points {
            let rij = p_i - pvec[j];
            let mut rij2 = (rij * rij).sum();
            if rij2 <= AC::H2 {
                rij2 = rij2.sqrt();
                let rij1 = 1.0 / (rij2+AC::EPS);
                let q = rij2 * AC::H1;
                let q1 = q - 2.0;
                let qq3 = q * q1 * q1 * q1;
                let wq = AC::AD * AC::FAC * qq3;

                wgx += wq * (rij.0 * rij1) * AC::H1;
                wgy += wq * (rij.1 * rij1) * AC::H1;
                wgz += wq * (rij.2 * rij1) * AC::H1;
            }
    }

    let u_i = uvec[iter];

    let bvt = -AC::BETA * AC::V * AC::DT;
    let zt  = -AC::ZETA * AC::DT;

    let dux = wgx*bvt + u_i.0*zt;
    let duy = wgy*bvt + u_i.1*zt;
    let duz = wgz*bvt + u_i.2*zt;

    let utmp_loop = u_i + point::Point(dux, duy, duz);
    let ptmp_loop = p_i + point::Point(dux * AC::DT, duy *  AC::DT, duz * AC::DT);

    UpdatePoints {utmp : utmp_loop,ptmp : ptmp_loop}
    
}