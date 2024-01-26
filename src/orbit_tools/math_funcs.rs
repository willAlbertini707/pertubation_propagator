// This module provides simple implementations of needed math functions


use std::f64::consts::PI;
use ndarray::{Array1, arr1};

pub fn deg2rad(degree: f64) -> f64 {
    // convert degrees to radians
    degree * PI / 180.0
}


pub fn cross(vect1: Array1<f64>, vect2: Array1<f64>) -> Array1<f64> {
    //implementation of cross product
    arr1(&[
        vect1[1]*vect2[2] - vect1[2] * vect2[1],
        vect1[2]*vect2[0] - vect1[0] * vect2[2],
        vect1[0]*vect2[1] - vect1[1] * vect2[0],
        ])
}



// ---------------------------------- tests -------------------------------
#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_cross_product() {

        let a = arr1(&[1.0, 2.0, 3.0]);
        let b = arr1(&[3.0, 2.0, 1.0]);

        let result = arr1(&[-4.0, 8.0, -4.0]);

        assert_eq!(result, cross(a,b));
    }
}