use point;

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
}

pub fn packstep_s(data: &mut point::Particles<'_>,
                  ncut: usize) {
    let n = data.pvec.len();

    for i in 0..(n - ncut) {
        packstep_single(data,i);
        data.uvec[i] = data.utmp[i];
        data.pvec[i] = data.ptmp[i];
    }
}


fn packstep_single(data: &mut point::Particles<'_>, iter: usize){
    let mut wgx = 0.0;
    let mut wgy = 0.0;
    let mut wgz = 0.0;
    let p_i = data.pvec[iter];

    for (j, p_j) in data.pvec.iter().enumerate() {
        if j != iter {
            let rij = p_i - *p_j;
            let mut rij2 = (rij * rij).sum();
            if rij2 <= AC::H2 {
                rij2 = rij2.sqrt();
                let rij1 = 1.0 / (rij2);
                let q = rij2 * AC::H1;
                let q1 = q - 2.0;
                let qq3 = q * q1 * q1 * q1;
                let wq = AC::AD * AC::FAC * qq3;

                wgx += wq * (rij.0 * rij1) * AC::H1;
                wgy += wq * (rij.1 * rij1) * AC::H1;
                wgz += wq * (rij.2 * rij1) * AC::H1;
            }
        }
    }

    let u_i = data.uvec[iter];

    let dux = (-AC::BETA * wgx * AC::V - AC::ZETA * u_i.0) * AC::DT;
    let duy = (-AC::BETA * wgy * AC::V - AC::ZETA * u_i.1) * AC::DT;
    let duz = (-AC::BETA * wgz * AC::V - AC::ZETA * u_i.2) * AC::DT;

    data.utmp[iter] = u_i + point::Point(dux, duy, duz);
    data.ptmp[iter] = p_i + point::Point(dux * AC::DT, duy, duz * AC::DT);
}