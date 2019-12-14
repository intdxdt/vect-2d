use super::*;
use geom_2d::{Point, Geometry, pt};
use crate::Vect;
use math_util::{round, PI};

const PREC: i32 = 8;

#[test]
fn test_zero_vector() {
//        let a = Point::new(0.88682, -1.06102);
//        let b = Point::new(3.5, 1.0);
//        let c = Point::new(-3.0, 1.0);
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
//        let b_pt = Point::new(3.5, 1.0);
//        let c_pt = Point::new(-3.0, 1.0);

    let a = [10., 150., 6.5];
    let e = [280., 280., 12.8];
    let v = Vect::new(a[0..2].into(), e[0..2].into());
    let pv = v.v;
    let nv = v.v.neg();
    let neg_a = Point::new(-a_pt.x, -a_pt.y);
    assert_eq!(nv, pv.kproduct(-1.0));
    assert_eq!(a_pt.neg(), neg_a);
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
//        let c_pt = Point::new(-3.0, 1.0);
    let u = Vect::new(Point::new(0., 0.), a_pt);
    let v = Vect::new(Point::new(0., 0.), b_pt);
    assert_eq!(round(u.project(&v), 5), 0.56121);
}

#[test]
fn test_vect() {
//        let a2 = pt![0.88682, -1.06102];
//        let b2 = pt![3.5, 1];
//        let c2 = pt![-3, 1];
//        let d2 = pt![-1.5, -3];

//        let m = 25.0;
//        let dir = (165.0f64).to_radians();
//        let cpt = Point::component(m, dir);
//        let va = Vect::new(Point::new(0., 0.), cpt);
//        let va_b = pt![-24.148145657226706, 6.470476127563026];

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

    assert_eq!(round(sed_v.magnitude(), PREC), 93.24400487);
    assert_eq!(round(sed_v2.magnitude(), PREC), 93.24400487);
}

#[test]
fn test_ext_vect() {
//        let a = Point::new(0.88682, -1.06102);
//        let b = Point::new(3.5, 1.0);
    let c = Point::new(-3.0, 1.0);

    let a2 = pt![0.88682, -1.06102];
    let b2 = pt![3.5, 1];
    let c2 = pt![-3, 1];
    let d2 = pt![-1.5, -3];

    let va = Vect::new(Point::new_origin(), a2);
    let vb = Vect::new(Point::new_origin(), b2);
    let vc = Vect::new(Point::new_origin(), c2);
    let vd = Vect::new(Point::new_origin(), d2);
    let vdb = Vect::new(d2, b2);

    assert_eq!(round(va.direction(), PREC),
               round(f64::to_radians(309.889497029295), PREC),
    );
    assert_eq!(round(vb.direction(), PREC),
               round(f64::to_radians(15.945395900922854), PREC),
    );
    assert_eq!(round(vc.direction(), PREC),
               round(f64::to_radians(161.565051177078), PREC),
    );
    assert_eq!(round(vd.direction(), PREC),
               round(f64::to_radians(243.43494882292202), PREC),
    );
    assert_eq!(va.a.x, 0.);
    assert_eq!(vc.a.x, vd.a.x);
    assert_eq!(round(vdb.magnitude(), 4),
               round(6.4031242374328485, 4),
    );
    assert_eq!(round(vdb.direction(), PREC),
               round(f64::to_radians(38.65980825409009), PREC),
    );
    let deflangle = 157.2855876468;
    let vo = vdb.extend_vect(3.64005494464026, f64::to_radians(180. + deflangle), true);
    let vo_defl = vdb.deflect_vector(3.64005494464026, f64::to_radians(-deflangle), true);
    // , "compare deflection and extending"
    assert_eq!(vo.b, vo_defl.b);
    // "vo by extending vdb by angle to origin"
    assert_eq!(round(vo.b.x, PREC), 0.0);
    // "vo by extending vdb by angle to origin"
    assert_eq!(round(vo.b.y, 4), round(0.0, PREC));
    let deflangle_b = 141.34019174590992;
    let inclangle_d = 71.89623696549336;
    // extend to c from end
    let vextc = vdb.extend_vect(6.5, f64::to_radians(180. + deflangle_b), true);
    ////extend to c from begining
    let vext_cfrom_d = vdb.extend_vect(4.272001872658765, f64::to_radians(inclangle_d), false);
    // deflect to c from begin
    let vdefl_cfrom_d = vdb.deflect_vector(4.272001872658765, f64::to_radians(180. - inclangle_d), false);
    // "comparing extend and deflect from begin point D"
    assert_eq!(vext_cfrom_d.b, vdefl_cfrom_d.b);
    // "vextc from b and from D : extending vdb by angle to c"
    assert_eq!(round(vext_cfrom_d.b.x, PREC), round(vextc.b.x, PREC));
    // "vextc from b and from D : extending vdb by angle to c"
    assert_eq!(round(vext_cfrom_d.b.y, PREC), round(vextc.b.y, PREC));
    // "vextc by extending vdb by angle to c"
    assert_eq!(round(vextc.b.x, PREC), c[0]);
    // "vextc by extending vdb by angle to c"
    assert_eq!(round(vextc.b.y, 4), c[1]);
    // "vextc with magnitudie extension from vdb c"
    assert_eq!(round(vextc.v.0.x, PREC), -vextc.magnitude());
    // "vextc horizontal vector test:  extension from vdb c"
    assert_eq!(round(vextc.v.0.y, PREC), 0.)
}

#[test]
fn test_vect_direction() {
    let m = 25.0;
    let dir = (165.0f64).to_radians();
    let cpt = Point::component(m, dir);
    let va = Vect::new(Point::new(0., 0.), cpt);
    let va_b = pt![-24.148145657226706, 6.470476127563026];

    //va endpoints equality: 0
    assert_eq!(round(va.b.x, PREC), round(va_b.x, PREC), );
    //va endpoints equality: 1
    assert_eq!(round(va.b.y, PREC), round(va_b.y, PREC), );
    assert!(va.magnitude().feq(25.));
    assert_eq!(f64::to_radians(165.), va.direction());
    assert_eq!(va.a.x, 0.0);
    assert_eq!(va.a.y, 0.0);

    //endpoint should be same as vector: 0
    assert_eq!(round(va.b.x, PREC), round(va.v.0.x, PREC));
    //endpoint should be same as vector: 1
    assert_eq!(round(va.b.y, PREC), round(va.v.0.y, PREC));
}


#[test]
fn test_vect_sideof() {
    let k = pt![-0.887, -1.6128];
    let u = pt![4.55309, 1.42996];

    let testpoints = vec![
        pt![2, 2], pt![0, 2], pt![0, -2], pt![2, -2], pt![0, 0], pt![2, 0], u, k,
    ];
    let v = Vect::new(k, u);
    let (mut left, mut right, mut on) = (Side::new(), Side::new(), Side::new());
    left.as_left();
    right.as_right();
    on.as_on();

    let mut sides = vec![Side::new(); testpoints.len()];
    for i in 0..testpoints.len() {
        sides[i] = v.side_of(testpoints[i]);
    }
    let side_out = vec![left, left, right, right, left, right, on, on];

    for i in 0..side_out.len() {
        assert!(sides[i].is_same_side(&side_out[i]))
    }
}
