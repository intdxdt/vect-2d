use geom_2d::{Point, Coordinate};
use side_rel::Side;
use math_util::Feq;
use std::cmp::Ordering;

#[derive(Copy, Clone, Debug)]
pub struct Vector(Point);

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

    ///Add creates A new point by adding to other point
    fn add(&self, o: Vector) -> Vector {
        Vector(self.0.add(&o.0))
    }

    ///Add creates A new point by adding to other point
    fn sub(&self, o: Vector) -> Vector {
        Vector(self.0.sub(&o.0))
    }

    ///Creates a new point by multiplying point by A scaler  k
    fn kproduct(&self, k: f64) -> Vector {
        Vector(self.0.kproduct(k))
    }

    ///Negate vector
    fn neg(&self) -> Vector {
        return self.kproduct(-1.0);
    }

    ///Computes vector magnitude of pt as vector: x , y as components
    fn magnitude(&self) -> f64 {
        self.0.magnitude()
    }

    ///Computes the square vector magnitude of pt as vector: x , y as components
    ///This has A potential overflow problem based on coordinates of pt x^2 + y^2
    fn square_magnitude(&self) -> f64 {
        self.0.square_magnitude()
    }

    ///Dot Product of two points as vectors
    fn dot_product(&self, o: Vector) -> f64 {
        self.0.dot_product(&o.0)
    }

    ///Unit vector of point
    fn unit_vector(&self) -> Vector {
        Vector(self.0.unit_vector())
    }

    ///project vector u on V
    fn project(&self, v: Vector) -> f64 {
        self.0.project(v.0)
    }

    ///Dir computes direction in radians - counter clockwise from x-axis.
    fn direction(&self) -> f64 {
        self.0.direction()
    }

    ///Reversed direction of vector direction
    fn reverse_direction(&self) -> f64 {
        geom_2d::reverse_direction(self.0.direction())
    }

    ///Computes the deflection angle from vector V to u
    fn deflection_angle(&self, u: Vector) -> f64 {
        geom_2d::deflection_angle(self.0.direction(), u.0.direction())
    }
}

impl Eq for Vector {}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.0.equals(&other.0)
    }
}

impl PartialOrd for Vector {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
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
        return self.v.0.magnitude();
    }

    ///Computes the direction of Vector
    pub fn direction(&self) -> f64 {
        return self.v.0.direction();
    }

    ///Reversed direction of vector direction
    pub fn reverse_direction(&self) -> f64 {
        return geom_2d::reverse_direction(self.v.0.direction());
    }

    ///Computes the deflection angle from vector V to u
    pub fn deflection_angle(&self, u: &Vect) -> f64 {
        return geom_2d::deflection_angle(self.direction(), u.direction());
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
mod tests {
    use super::*;
    use geom_2d::{Point, Geometry, Points, pt};
    use crate::Vect;
    use math_util::{round, PI};

    const prec: i32 = 8;

    #[test]
    fn test_zero_vector() {
        let A = Point::new(0.88682, -1.06102);
        let B = Point::new(3.5, 1.0);
        let C = Point::new(-3.0, 1.0);
        let v = Vect::new(Point::new_origin(), Point::new_origin());
        assert_eq!(v.a, Point::new(0., 0.));
        assert_eq!(v.b, Point::new(0., 0.));
        assert_eq!(v.v, Vector(Point::new(0., 0.)));
        let mdatbt = [v.magnitude(), v.direction(), v.at, v.bt];
        for o in mdatbt.iter() {
            assert_eq!(*o, 0.0);
        }
    }


    #[test]
    fn test_vector_component() {
        let v = Vect::new(Point::new_origin(), Point::new(3., 4.));
        assert_eq!(v.a, Point::new_origin());
        assert_eq!(v.b, Point::new(3., 4.));
        assert_eq!(v.v, Vector(Point::new(3., 4.)));
        assert_eq!(v.magnitude(), 5.0);
    }

    #[test]
    fn test_vector_negation() {
        let a_pt = Point::new(0.88682, -1.06102);
        let b_pt = Point::new(3.5, 1.0);
        let c_pt = Point::new(-3.0, 1.0);

        let a = [10., 150., 6.5];
        let e = [280., 280., 12.8];
        let v = Vect::new(a[0..2].into(), e[0..2].into());
        let pv = v.v;
        let nv = v.v.neg();
        let mut negA = Point::new(-a_pt.x, -a_pt.y);
        assert_eq!(nv, pv.kproduct(-1.0));
        assert_eq!(a_pt.neg(), negA);
    }

    #[test]
    fn test_vect_construct() {
        let a = [10., 150., 6.5];
        let e = [280., 280., 12.8];
        let i = [185., 155., 8.6];
        let v = Vect::new_with_t(a[..2].into(), e[..2].into(), a[2], e[2]);
        let vo = Vect::new_with_t(a[..2].into(), e[..2].into(), a[2], e[2]);
        let vi = Vect::new_with_t(i[..2].into(), e[..2].into(), i[2], e[2]);
        let m = 5.0;
        let dir = 53.13010235415598f64.to_radians();
        let c_pt = Point::component(m, dir);
        let vk = Vect::new(Point::new_origin(), c_pt);

        assert_eq!(vk.a, Point::new_origin());
        assert_eq!(round(vk.b.x, 8), 3.0);
        assert_eq!(round(vk.b.y, 8), 4.0);

        assert_eq!(v.a, vo.a);
        assert_eq!(v.a, vo.a);
        assert_eq!(v.b, vo.b);
        assert_eq!(v.b, vo.b);

        assert_eq!(v.magnitude(), vo.magnitude());
        assert_eq!(v.magnitude(), vo.magnitude());
        assert_eq!(v.direction(), vo.direction());
        assert_eq!(v.direction(), vo.direction());

        assert_eq!(v.a.as_array(), a[0..2]);
        assert_eq!(v.b.as_array(), e[0..2]);
        assert_eq!(vi.a.as_array(), [i[0], i[1]]);
        assert_eq!(vi.b, v.b);

        assert_eq!(v.at, a[2]);
        assert_eq!(v.bt, e[2]);
        assert_eq!(v.dt(), e[2] - a[2]);

        let aa = Point::new(a[0], a[1]);
        let ee = Point::new(e[0], e[1]);
        let d = ee.distance(&aa);
        assert_eq!(v.magnitude(), d);
    }

    #[test]
    fn test_direction() {
        let v = Vect::new(Point::new_origin(), Point::new(-1., 0.));
        assert_eq!(v.direction(), PI);
        assert_eq!(Vector([1., 1.].into()).direction(), std::f64::consts::FRAC_PI_4);
        assert_eq!(Vector([-1., 0.].into()).direction(), PI);
        assert_eq!(Vector([1., 3f64.sqrt()].into()).direction(), 60f64.to_radians());
        assert_eq!(Vector([0, -1].into()).direction(), 270f64.to_radians())
    }

    #[test]
    fn test_rev_direction() {
        let v = Vect::new(Point::new_origin(), Point::new(-1., 0.));
        assert_eq!(v.direction(), PI);
        assert_eq!(v.reverse_direction(), 0.0)
    }

    #[test]
    fn test_deflection_angle() {
        let ln0 = [pt![0, 0], pt![20, 30]];
        let ln1 = [pt![20, 30], pt![40, 15]];
        let v0 = Vect::new(ln0[0], ln0[1]);
        let v1 = Vect::new(ln1[0], ln1[1]);

        assert_eq!(
            round(v0.deflection_angle(&v1), 10),
            round(93.17983011986422f64.to_radians(), 10)
        );
        assert_eq!(round(v0.deflection_angle(&v0), 10), 0f64.to_radians());

        let ln1 = [pt![20, 30], pt![20, 60]];
        let v1 = Vect::new(ln1[0], ln1[1]);

        assert_eq!(
            round(v0.deflection_angle(&v1), 10),
            round(-33.690067525979806f64.to_radians(), 10)
        );
    }

    #[test]
    fn test_projection() {
        let a_pt = Point::new(0.88682, -1.06102);
        let b_pt = Point::new(3.5, 1.0);
        let c_pt = Point::new(-3.0, 1.0);
        let u = Vect::new(Point::new(0., 0.), a_pt);
        let v = Vect::new(Point::new(0., 0.), b_pt);
        assert_eq!(round(u.project(&v), 5), 0.56121);
    }

    #[test]
    fn test_vect() {
        let A2 = pt![0.88682, -1.06102];
        let B2 = pt![3.5, 1];
        let C2 = pt![-3, 1];
        let D2 = pt![-1.5, -3];

        let m = 25.0;
        let dir = (165.0f64).to_radians();
        let cpt = Point::component(m, dir);
        let va = Vect::new(Point::new(0., 0.), cpt);
        let va_b = pt![-24.148145657226706, 6.470476127563026];

        //should test distance vector
        let a = pt![16.82295, 10.44635];
        let b = pt![28.99656, 15.76452];
        let on_ab = pt![25.32, 14.16];

        let tpoints = vec![
            pt![30., 0.], pt![15.78786, 25.26468], pt![-2.61504, -3.09018], pt![28.85125, 27.81773],
            a, b, on_ab,
        ];

        let t_dists = vec![14.85, 13.99, 23.69, 12.05, 0.00, 0.00, 0.00];
        let tvect = Vect::new(a, b);
        let mut dists = vec![0f64; tpoints.len()];
        for i in 0..tpoints.len() {
            let tp = tpoints[i];
            dists[i] = tvect.distance_to_point(&tp);
        }

        for i in 0..tpoints.len() {
            assert_eq!(round(dists[i], 2), round(t_dists[i], 2))
        }
    }

    #[test]
    fn test_sed_vect() {
        //should test side sed vector to point at time T
        let a = [10., 150., 6.5];
        let e = [280., 280., 12.8];
        let i = [185., 155., 8.6];
        let ai: Point = [i[0], i[1]].into();
        let v = Vect::new_with_t(
            a[..2].into(), e[..2].into(), a[2], e[2],
        );

        let sed_v = v.sed_vector(ai, i[2]);
        let sed_v2 = v.sed_vector(ai, i[2]);

        assert_eq!(round(sed_v.magnitude(), prec), 93.24400487);
        assert_eq!(round(sed_v2.magnitude(), prec), 93.24400487);
    }

    #[test]
    fn test_ext_vect() {
        let A = Point::new(0.88682, -1.06102);
        let B = Point::new(3.5, 1.0);
        let C = Point::new(-3.0, 1.0);

        let A2 = pt![0.88682, -1.06102];
        let B2 = pt![3.5, 1];
        let C2 = pt![-3, 1];
        let D2 = pt![-1.5, -3];

        let va = Vect::new(Point::new_origin(), A2);
        let vb = Vect::new(Point::new_origin(), B2);
        let vc = Vect::new(Point::new_origin(), C2);
        let vd = Vect::new(Point::new_origin(), D2);
        let vdb = Vect::new(D2, B2);

        assert_eq!(round(va.direction(), prec),
                   round(f64::to_radians(309.889497029295), prec),
        );
        assert_eq!(round(vb.direction(), prec),
                   round(f64::to_radians(15.945395900922854), prec),
        );
        assert_eq!(round(vc.direction(), prec),
                   round(f64::to_radians(161.565051177078), prec),
        );
        assert_eq!(round(vd.direction(), prec),
                   round(f64::to_radians(243.43494882292202), prec),
        );
        assert_eq!(va.a.x, 0.);
        assert_eq!(vc.a.x, vd.a.x);
        assert_eq!(round(vdb.magnitude(), 4),
                   round(6.4031242374328485, 4),
        );
        assert_eq!(round(vdb.direction(), prec),
                   round(f64::to_radians(38.65980825409009), prec),
        );
        let deflangle = 157.2855876468;
        let vo = vdb.extend_vect(3.64005494464026, f64::to_radians(180. + deflangle), true);
        let vo_defl = vdb.deflect_vector(3.64005494464026, f64::to_radians(-deflangle), true);
        // , "compare deflection and extending"
        assert_eq!(vo.b, vo_defl.b);
        // "vo by extending vdb by angle to origin"
        assert_eq!(round(vo.b.x, prec), 0.0);
        // "vo by extending vdb by angle to origin"
        assert_eq!(round(vo.b.y, 4), round(0.0, prec));
        let deflangleB = 141.34019174590992;
        let inclangleD = 71.89623696549336;
        // extend to c from end
        let vextc = vdb.extend_vect(6.5, f64::to_radians(180. + deflangleB), true);
        ////extend to c from begining
        let vextCFromD = vdb.extend_vect(4.272001872658765, f64::to_radians(inclangleD), false);
        // deflect to c from begin
        let vdeflCFromD = vdb.deflect_vector(4.272001872658765, f64::to_radians(180. - inclangleD), false);
        // "comparing extend and deflect from begin point D"
        assert_eq!(vextCFromD.b, vdeflCFromD.b);
        // "vextc from B and from D : extending vdb by angle to C"
        assert_eq!(round(vextCFromD.b.x, prec), round(vextc.b.x, prec));
        // "vextc from B and from D : extending vdb by angle to C"
        assert_eq!(round(vextCFromD.b.y, prec), round(vextc.b.y, prec));
        // "vextc by extending vdb by angle to C"
        assert_eq!(round(vextc.b.x, prec), C[0]);
        // "vextc by extending vdb by angle to C"
        assert_eq!(round(vextc.b.y, 4), C[1]);
        // "vextc with magnitudie extension from vdb C"
        assert_eq!(round(vextc.v.0.x, prec), -vextc.magnitude());
        // "vextc horizontal vector test:  extension from vdb C"
        assert_eq!(round(vextc.v.0.y, prec), 0.)
    }
}
