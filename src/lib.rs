use geom_2d::{Point, Coordinate};
use side_rel::Side;
use math_util::Feq;

struct Vector(Point);

impl Vector {
    pub fn new(a: Point, b: Point) -> Vector {
        Vector(b.sub(&a))
    }
    pub fn from_pt(pt: Point) -> Vector {
        Vector(pt)
    }
    pub fn from_xy(x: f64, y: f64) -> Vector {
        Vector(Point::new(x, y))
    }
}

struct Vect {
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

    ///Magnitude of Vector
    fn magnitude(&self) -> f64 {
        return self.v.0.magnitude();
    }

    ///Computes the Direction of Vector
    fn direction(&self) -> f64 {
        return self.v.0.direction();
    }

    ///Reversed direction of vector direction
    pub fn reverse_direction(&self) -> f64 {
        return geom_2d::reverse_direction(self.v.0.direction());
    }


    ///Computes the deflection angle from vector V to u
    fn deflection_angle(&self, u: &Vect) -> f64 {
        return geom_2d::deflection_angle(self.direction(), u.direction());
    }

    ///computes the change in time between bt and at
    fn dt(&self) -> f64 {
        return (self.bt - self.at).abs();
    }


    ///Computes the relation of a point to vector
    fn side_of(&self, pnt: Point) -> Side {
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
    fn sed_vector(&self, pnt: Point, t: f64) -> Vect {
        let m = (self.magnitude() / self.dt()) * (t - self.at);
        //var vb = v.extend_vect(m, 0.0, false)
        let c_pt = self.v.0.extend(m, 0., false);
        let a = self.a.add(&c_pt);
        return Vect::new(a, pnt);
    }

    ///Extends vector from the from end or from begin of vector
    fn extend_vect(&self, magnitude: f64, angle: f64, from_end: bool) -> Vect {
        let c_pt = self.v.0.extend(magnitude, angle, from_end);
        let cv = Vector::from_pt(c_pt);
        let a = if from_end { self.b } else { self.a };
        Vect { a, b: a.add(&cv.0), v: cv, at: 0., bt: 0. }
    }

    ///Computes vector deflection given deflection angle and
    /// side of vector to deflect from (from_end)
    fn deflect_vector(&self, magnitude: f64, defl_angle: f64, from_end: bool) -> Vect {
        let c_pt = self.v.0.deflect(magnitude, defl_angle, from_end);
        let cv = Vector::from_pt(c_pt);
        let a = if from_end { self.b } else { self.a };
        Vect { a, b: a.add(&cv.0), v: cv, at: 0., bt: 0. }
    }

    ///Computes distance from A point to Vect
    fn distance_to_point(&self, pnt: &Point) -> f64 {
        geom_2d::segment::distance_to_point(&self.a, &self.b, pnt)
    }

    ///Project vector u on V
    fn project(&self, onv: &Vect) -> f64 {
        self.v.0.project(onv.v.0)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
