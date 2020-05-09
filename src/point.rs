use std::fmt;
use std::ops;

#[derive(Copy, Clone)]
pub struct Point(pub f32, pub f32, pub f32);

pub struct Particles{
    pub pvec: Vec<Point>,
    pub uvec: Vec<Point>,
}

// Implement method for print (println!)
impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, {})", self.0, self.1, self.2)
    }
}

// Overload operators for "Point"
impl ops::Add for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point(self.0 + other.0, self.1 + other.1, self.2 + other.2)
    }
}

impl ops::Sub for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point(self.0 - other.0, self.1 - other.1, self.2 - other.2)
    }
}

impl ops::Div for Point {
    type Output = Point;

    fn div(self, other: Point) -> Point {
        Point(self.0 / other.0, self.1 / other.1, self.2 / other.2)
    }
}

impl ops::Mul for Point {
    type Output = Point;

    fn mul(self, other: Point) -> Point {
        Point(self.0 * other.0, self.1 * other.1, self.2 * other.2)
    }
}

// Implement methods on struct
impl Point {
    pub fn new() -> Self {
        Point(0.0, 0.0, 0.0)
    }

    // Sum
    pub fn sum(&self) -> f32 {
        self.0 + self.1 + self.2
    }
}
