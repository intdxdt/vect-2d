use geom_2d::{Point, Coordinate};
use side_rel::Side;
use math_util::Feq;

struct Vector(Point);

impl Vector {
    pub fn new(a: Point, b: Point) -> Vector {
        Vector(b.sub(&a))
    }
    pub fn new_xy(x: f64, y: f64) -> Vector {
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
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
