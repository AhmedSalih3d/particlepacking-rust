use point;
use std::fs::File;
use std::io::{BufRead, BufReader};

// Read data file
// https://stackoverflow.com/questions/29580979/how-to-read-a-matrix-from-a-txt-file-in-rust
pub fn readpoints(filename: &str) -> Vec<point::Point> {
    let f = BufReader::new(File::open(filename).unwrap());

    // Basically skip first line
    //f.read_line(&mut s).unwrap();

    let arr: Vec<Vec<f32>> = f
        .lines()
        .map(|l| {
            l.unwrap()
                .split(char::is_whitespace)
                .map(|number| number.parse().unwrap())
                .collect()
        })
        .collect();

    let n = arr.len();

    let mut pvec = vec![point::Point::new(); n];
    for i in 0..n {
        pvec[i] = point::Point(arr[i][0], arr[i][1], arr[i][2]);
    }
    pvec
}
