use crate::Vector;
use math_util::{round, Feq};
use geom_2d::{Point, pt};

#[test]
fn test_vector() {
    //should test unit vector
    let a = Vector::from_xy(237., 289.);
    let b = Vector::from_xy(462., 374.);
    let ab = a.add(b);
    assert_eq!(ab[0], 237. + 462.);
    assert_eq!(ab[1], 289. + 374.);

    //let v0 = Vector::from_xy(0., 0.);
    let v = Vector::from_xy(-3., 2.);
    let mut unit_v = v.unit_vector();
    unit_v.0.x = round(unit_v.0.x, 6);
    unit_v.0.y = round(unit_v.0.y, 6);
    assert_eq!(unit_v.0, Point::new(-0.83205, 0.5547));

    let v3 = Vector::new(pt![0., 0.], pt![3., 4.]);
    assert_eq!(v3.magnitude(), 5.0);
    let v4 = v3.sub(Vector::from_xy(2., 3.));
    assert_eq!(v4.square_magnitude(), 2.0);
    assert_eq!(v4.magnitude(), math_util::SQRT_2);


    //Vector - unit & Project
    let u = Vector::from_xy(0.88682, -1.06102);
    let v = Vector::from_xy(3.5, 1.0);
    assert_eq!(round(u.project(v), 5), 0.56121);
    let z = Vector::from_xy(0., 0.);
    let zv = z.unit_vector();
    assert!(zv[0].feq(0.));
    assert!(zv[1].feq(0.));

    //should test dot product
    let dot_prod = Vector::from_xy(1.2, -4.2).dot_product(
        Vector::from_xy(1.2, -4.2),
    );
    assert_eq!(19.08, round(dot_prod, 8))
}
