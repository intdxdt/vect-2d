use geom_2d::{Point, Coordinate};
use side_rel::Side;
use math_util::Feq;
use std::ops::Index;

#[derive(Copy, Clone, Debug)]
pub struct Vector(Point);

impl Vector {
    ///New Vector constructor
    pub fn new(a: Point, b: Point) -> Vector {
        Vector(b.sub(&a))
    }

    ///Vector from point
    pub fn from_pt(pt: Point) -> Vector {
        Vector(pt)
    }

    ///Vector from x, y components
    pub fn from_xy(x: f64, y: f64) -> Vector {
        Vector(Point::new(x, y))
    }

    ///Add creates A new point by adding to other point
    pub fn add(&self, o: Vector) -> Vector {
        Vector(self.0.add(&o.0))
    }

    ///Add creates A new point by adding to other point
    pub fn sub(&self, o: Vector) -> Vector {
        Vector(self.0.sub(&o.0))
    }

    ///Creates a new point by multiplying point by A scaler  k
    pub fn kproduct(&self, k: f64) -> Vector {
        Vector(self.0.kproduct(k))
    }

    ///Negate vector
    pub fn neg(&self) -> Vector {
        return self.kproduct(-1.0);
    }

    ///Computes vector magnitude of pt as vector: x , y as components
    pub fn magnitude(&self) -> f64 {
        self.0.magnitude()
    }

    ///Computes the square vector magnitude of pt as vector: x , y as components
    ///This has A potential overflow problem based on coordinates of pt x^2 + y^2
    pub fn square_magnitude(&self) -> f64 {
        self.0.square_magnitude()
    }

    ///Dot Product of two points as vectors
    pub fn dot_product(&self, o: Vector) -> f64 {
        self.0.dot_product(&o.0)
    }

    ///Unit vector of point
    pub fn unit_vector(&self) -> Vector {
        Vector(self.0.unit_vector())
    }

    ///project vector u on V
    pub fn project(&self, v: Vector) -> f64 {
        self.0.project(v.0)
    }

    ///Dir computes direction in radians - counter clockwise from x-axis.
    pub fn direction(&self) -> f64 {
        self.0.direction()
    }

    ///Reversed direction of vector direction
    pub fn reverse_direction(&self) -> f64 {
        geom_2d::reverse_direction(self.0.direction())
    }

    ///Computes the deflection angle from vector V to u
    pub fn deflection_angle(&self, u: Vector) -> f64 {
        geom_2d::deflection_angle(self.0.direction(), u.0.direction())
    }
}

impl Index<usize> for Vector {
    type Output = f64;
    fn index(&self, i: usize) -> &Self::Output {
        match i {
            0 => &self.0.x,
            1 => &self.0.y,
            _ => unreachable!(),
        }
    }
}

impl Eq for Vector {}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.0.equals(&other.0)
    }
}


pub struct Vect {
    pub a: Point,
    pub b: Point,
    pub v: Vector,
    pub at: f64,
    pub bt: f64,
}

impl Vect {
    ///New vect given a and b
    pub fn new(a: Point, b: Point) -> Vect {
        Vect { a, b, v: Vector::new(a, b), at: 0., bt: 0. }
    }
    ///New vect given a,b,at & bt
    pub fn new_with_t(a: Point, b: Point, at: f64, bt: f64) -> Vect {
        Vect { a, b, v: Vector::new(a, b), at, bt }
    }

    ///magnitude of Vector
    pub fn magnitude(&self) -> f64 {
        return self.v.magnitude();
    }

    ///Computes the direction of Vector
    pub fn direction(&self) -> f64 {
        return self.v.direction();
    }

    ///Reversed direction of vector direction
    pub fn reverse_direction(&self) -> f64 {
        return self.v.reverse_direction();
    }

    ///Computes the deflection angle from vector V to u
    pub fn deflection_angle(&self, u: &Vect) -> f64 {
        return self.v.deflection_angle(u.v);
    }

    ///computes the change in time between bt and at
    pub fn dt(&self) -> f64 {
        return (self.bt - self.at).abs();
    }

    ///Computes the relation of a point to vector
    pub fn side_of(&self, pnt: Point) -> Side {
        let ccw = pnt.orientation2d(&self.a, &self.b);
        let mut s = Side::new();
        s.as_left();
        if ccw.feq(0.) {
            s.as_on();
        } else if ccw > 0. {
            s.as_right();
        }
        s
    }

    ///Computes the Synchronized Euclidean Distance - Vector
    pub fn sed_vector(&self, pnt: Point, t: f64) -> Vect {
        let m = (self.magnitude() / self.dt()) * (t - self.at);
        //var vb = v.extend_vect(m, 0.0, false)
        let c_pt = self.v.0.extend(m, 0., false);
        let a = self.a.add(&c_pt);
        return Vect::new(a, pnt);
    }

    ///Extends vector from the from end or from begin of vector
    pub fn extend_vect(&self, magnitude: f64, angle: f64, from_end: bool) -> Vect {
        let c_pt = self.v.0.extend(magnitude, angle, from_end);
        let cv = Vector::from_pt(c_pt);
        let a = if from_end { self.b } else { self.a };
        Vect { a, b: a.add(&cv.0), v: cv, at: 0., bt: 0. }
    }

    ///Computes vector deflection given deflection angle and
    /// side of vector to deflect from (from_end)
    pub fn deflect_vector(&self, magnitude: f64, defl_angle: f64, from_end: bool) -> Vect {
        let c_pt = self.v.0.deflect(magnitude, defl_angle, from_end);
        let cv = Vector::from_pt(c_pt);
        let a = if from_end { self.b } else { self.a };
        Vect { a, b: a.add(&cv.0), v: cv, at: 0., bt: 0. }
    }

    ///Computes distance from A point to Vect
    pub fn distance_to_point(&self, pnt: &Point) -> f64 {
        geom_2d::segment::distance_to_point(&self.a, &self.b, pnt)
    }

    ///project vector u on V
    pub fn project(&self, onv: &Vect) -> f64 {
        self.v.0.project(onv.v.0)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod tests_vector;