mod algorithm;
mod dataio;
mod plotting;
mod point;

extern crate charts;
extern crate rayon;
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

fn rescale_data2d(pvec: &mut Vec<point::Point>){
    let minx = pvec.iter().fold(0.0f32, |min, &val| if val.0 < min{ val.0 } else{ min });
    let maxx = pvec.iter().fold(0.0f32, |max, &val| if val.0 > max{ val.0 } else{ max });
    let minz = pvec.iter().fold(0.0f32, |min, &val| if val.2 < min{ val.2 } else{ min });
    let maxz = pvec.iter().fold(0.0f32, |max, &val| if val.2 > max{ val.2 } else{ max });
    //pvec.iter().fold(0.0f32, |min_val, &val| val.0.min(min_val));
    for e in pvec {
        *e = *e + point::Point(-minx,0.0,-minz);
        *e = *e / point::Point( maxx,1.0,maxz);
    }

}

fn main() {
    rayon::ThreadPoolBuilder::new().num_threads(8).build_global().unwrap();
    const FILENAME:&str = "RustCircleUniformGridSQUARE.txt";
    let mut pini = dataio::readpoints(FILENAME);
    let uini = genvelocity(&pini);
    rescale_data2d(&mut pini);

    let mut data = point::Particles { pvec: pini, 
                                      uvec: uini,  
                                };

                                

    const NCUT: usize = 800;
    let now = Instant::now();
    for _ in 0..100 {
        algorithm::packstep_s(&mut data, NCUT);
    }
    println!("{}", now.elapsed().as_millis());

    const NAME:&str = "PointsYay.svg";
    // How to transfer ownership correctly..
    plotting::plot_points(&data.pvec, NAME);
}
