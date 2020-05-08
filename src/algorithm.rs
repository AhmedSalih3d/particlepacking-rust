use point;

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

pub fn packstep_s(data: &mut point::Particles<'_>,
                  ncut: usize) {
    let n = data.pvec.len();
    let pi_vec = data.pvec[0..(n-ncut)].to_vec();

    let pim = &*data.pvec.clone();
    let uim = &*data.uvec.clone();

    for (i,_) in pi_vec.iter().enumerate() {
        let up_points = packstep_single(&pim,&uim,i);
        data.ptmp[i] = up_points.ptmp;
        data.utmp[i] = up_points.utmp;
    }


    data.pvec[..(n - ncut)].clone_from_slice(&data.ptmp[..(n - ncut)]);
    data.uvec[..(n - ncut)].clone_from_slice(&data.utmp[..(n - ncut)]);
}

fn packstep_single(pvec: &[point::Point], uvec: & &[point::Point], iter: usize) -> UpdatePoints{
    let mut wgx = 0.0;
    let mut wgy = 0.0;
    let mut wgz = 0.0;
    let p_i = pvec[iter];

    for p_j in pvec.iter() {
            let rij = p_i - *p_j;
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