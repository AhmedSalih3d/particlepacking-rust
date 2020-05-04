mod algorithm;
mod dataio;
mod plotting;
mod point;

extern crate charts;
use std::time::Instant;

fn genvelocity(pvec: &[point::Point]) -> Vec<point::Point> {
    let a0 = 50.0;
    let n = pvec.len();
    let mut uvec = vec![point::Point::new(); n];
    for i in 0..n {
        uvec[i] = pvec[i] * point::Point(a0, 0.0, -a0);
    }
    uvec
}

fn main() {
    const FILENAME:&str = "RustCircleUniformGridSQUARE.txt";
    let mut pini = dataio::readpoints(FILENAME);
    let mut uini = genvelocity(&pini);
    let mut ptmp = pini.clone();
    let mut utmp = uini.clone();

    let mut data = point::Particles { pvec: &mut pini, 
                                      uvec: &mut uini,  
                                      ptmp: &mut ptmp,
                                      utmp: &mut utmp
                                };

    const NCUT: usize = 800;
    let now = Instant::now();
    for _ in 0..100 {
        algorithm::packstep_s(&mut data, NCUT);
    }
    println!("{}", now.elapsed().as_millis());

    const NAME:&str = "Points.svg";
    // How to transfer ownership correctly..
    plotting::plot_points(&pini, NAME);
}
