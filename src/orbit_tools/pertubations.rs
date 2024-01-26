// external imports
use ndarray::{Array1, arr1};

// internal imports
use super::propagator::OrbitSystem;
use super::math_funcs::{cross, deg2rad};


// ----------------------------------- ATMOSHPERIC DRAG ---------------------------------------------
fn exp_rho(altitude: f64) -> f64 {

	let mut z = altitude;
	// exponential atmosphere model implementation from Orbital Mechanics
    // for Engineering Students (Howard D. Curtis)

	let h = [ 0.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0,
					150.0, 180.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0];

	
	let r = [1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5, 
	 1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9, 
	 2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12,
	 1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15];

	let H = [7.310, 6.427, 6.546, 7.360, 8.342, 7.583, 6.661, 
	 5.927, 5.533, 5.703, 6.782, 9.973, 13.243, 16.322, 
	 21.652, 27.974, 34.934, 43.342, 49.755, 54.513, 58.019, 
	 60.980, 65.654, 76.377, 100.587, 147.203, 208.020];

	let mut i: usize = 0;

	if z > 1000.0 {
        z = 1000.0;
    } else if z < 0.0{
        z = 0.0;
    }

	for j in 0..27 {

		if z >= h[j] && z < h[j + 1] {
            i = j;
        }

		if z == 1000.0 {
            i = 26;
        }
    }
	
	r[i]* (-(z - h[i])/H[i]).exp()
}

pub fn drag_acceleration(r: Array1<f64>, r_mag: f64, v: Array1<f64>, vals: &OrbitSystem) -> Array1<f64> {
    // drag acceleration model

    let v_rel = v - cross(arr1(&[0.0, 0.0, 72.9211e-6]), r);
    let v_rel_norm = v_rel.dot(&v_rel).sqrt();

    return -0.5 * vals.cd * 1000.0 * vals.area * exp_rho(r_mag-6378.0) / vals.mass * v_rel * v_rel_norm;
}

// ----------------------------------- N-body and SRP ---------------------------------------------
pub fn solar_position(jd: f64) -> Array1<f64> {
    //
    // This function alculates the geocentric equatorial position vector
    // of the sun, given the julian date.
    //
    // User M-functions required: None
    // -------------------------------------------------------------------------
    //...Astronomical unit (km):
    let AU = 149597870.691;
    //...Julian days since J2000:
    let n = jd - 2451545.0;
    //...Julian centuries since J2000:
    let cy = n/36525.0;
    //...Mean anomaly (deg):
    let mut M = 357.528 + 0.9856003*n;
    M = M.rem_euclid(360.0);
    //...Mean longitude (deg):
    let mut L = 280.460 + 0.98564736*n;
    L = L.rem_euclid(360.0);
    //...Apparent ecliptic longitude (deg):
    let mut lamda = L + 1.915 * deg2rad(M).sin() + 0.020*deg2rad(2.0 * M).sin();
    lamda = lamda.rem_euclid(360.0);
    //...Obliquity of the ecliptic (deg):
    let eps = 23.439 - 0.0000004*n;
    //...Unit vector from earth to sun:
    let u = arr1(&[deg2rad(lamda).cos(),
                  deg2rad(lamda).sin()*deg2rad(eps).cos(), 
                  deg2rad(lamda).sin()*deg2rad(eps).sin()]);
    //...Distance from earth to sun (km):
    let rS = (1.00014 - 0.01671*deg2rad(M).cos() - 0.000140*deg2rad(2.0*M).cos()) * AU;


    //...Geocentric position vector (km):
    rS*u
}

pub fn three_body_solar(r: Array1<f64>, r_sun: Array1<f64>, r_sun_mag: f64) -> Array1<f64> {
    // gravtational constant of sun
    let mu_sun = 132712e6;

    // spacecraft relative to sun
    let r_sc_sun = r_sun.clone() - r;
    let r_sc_sun_mag = r_sc_sun.dot(&r_sc_sun).sqrt();

    // calculate n body acceleration
    let n_accel = mu_sun * (r_sc_sun/r_sc_sun_mag.powf(3.0) - r_sun/r_sun_mag.powf(3.0));
    
    n_accel
}

pub fn srp(r: Array1<f64>, r_mag: f64, r_sun: Array1<f64>, r_sun_mag: f64, system_vals: &OrbitSystem) -> Array1<f64> {
    // srp model
    let p_sr=4.57e-6;

    let theta_A = (system_vals.body_radius / r_sun_mag).acos();
    let theta_B = (system_vals.body_radius/r_mag).acos();
    let theta =  (r_sun.dot(&r) / (r_mag * r_sun_mag)).acos();

    let mut F = 1.0;

    if theta_A + theta_B < theta {
        F=0.0;
    }

    -p_sr * system_vals.cr * system_vals.area * r / system_vals.mass / r_mag * F / 1000.0
}

pub fn j2_j6(r: Array1<f64>, r_mag: f64, system_vals: &OrbitSystem) -> Array1<f64> {
    // calculate j2 - j6 effects
    let vals = system_vals;
    let ri = r[0];
    let rj = r[1];
    let rk = r[2];

    // define constants
    let J2 = 1.08262668355e-3;
    let J3 = -2.53265648533e-6;
    let J4 = -1.61962215937e-6;
    let J5 = -2.27296082869e-7;
    let J6 = 5.40681239107e-7;

    // j2
    let mut front_val = -3.0*J2*vals.mu*vals.body_radius.powf(2.0) / (2.0*r_mag.powf(5.0));
    let mut  back_val = 5.0 * rk.powf(2.0) / r_mag.powf(2.0);
    let j2_x: f64 = front_val * ri * (1.0 - back_val);
    let j2_y: f64 = front_val * rj * (1.0 - back_val);
    let j2_z: f64 = front_val * rk * (3.0 - back_val);

    // j3
    front_val = -5.0 * J3 * vals.mu * vals.body_radius.powf(3.0) / (2.0*r_mag.powf(7.0));
    back_val = 3.0 * rk - 7.0 * rk.powf(3.0) / r_mag.powf(2.0);
    let j3_x: f64 = front_val * ri * back_val;
    let j3_y: f64 = front_val * rj * back_val;
    let j3_z: f64 = front_val  * (6.0*rk.powf(2.0) - 7.0*rk.powf(4.0)/r_mag.powf(2.0) - 3.0* r_mag.powf(2.0)/5.0);

    // j4
    front_val = 15.0 * J4 * vals.mu * vals.body_radius.powf(4.0) / (8.0 * r_mag.powf(7.0));
    back_val = 1.0 - 14.0 * rk.powf(2.0) / r_mag.powf(2.0) + 21.0 * rk.powf(4.0) / r_mag.powf(4.0);
    let j4_x: f64 = front_val * ri * back_val;
    let j4_y: f64 = front_val * rj * back_val;
    let j4_z: f64 = front_val * rk * (5.0 - 70.0 * rk.powf(2.0) / 3.0 / r_mag.powf(2.0) + 21.0 * rk.powf(4.0) / r_mag.powf(4.0));

    // j5
    front_val = 3.0*J5*vals.mu*vals.body_radius.powf(5.0) / (8.0*r_mag.powf(9.0));
    back_val = 35.0 - 210.0 * rk.powf(2.0)/r_mag.powf(2.0) + 231.0 *rk.powf(4.0) / r_mag.powf(4.0);
    let j5_x: f64 = front_val * ri * rk * back_val;
    let j5_y: f64 = front_val * rj * rk * back_val;
    let j5_z: f64 = front_val * rk * rk * (105.0 - 315.0*rk.powf(2.0)/r_mag.powf(2.0) + 231.0*rk.powf(4.0)/r_mag.powf(4.0))
                     - 15.0 * J5 * vals.mu * vals.body_radius.powf(5.0) / (8.0 * r_mag.powf(7.0));

    // j6
    front_val = -J6*vals.mu*vals.body_radius.powf(6.0) / (16.0*r_mag.powf(9.0));
    back_val = 35.0 - 945.0*rk.powf(2.0)/r_mag.powf(2.0) + 3465.0*rk.powf(4.0)/r_mag.powf(4.0) - 3003.0*rk.powf(6.0)/r_mag.powf(6.0);
    let j6_x: f64 = front_val * ri * back_val;
    let j6_y: f64 = front_val * rj * back_val;
    let j6_z: f64 = front_val * rk * (245.0 - 2205.0*rk.powf(2.0)/r_mag.powf(2.0) + 
                                    4851.0 * rk.powf(4.0) / r_mag.powf(4.0) - 3003.0 * rk.powf(6.0) / r_mag.powf(6.0));

    // return sum
    arr1(&[
        j2_x + j3_x + j4_x + j5_x + j6_x,
        j2_y + j3_y + j4_y + j5_y + j6_y,
        j2_z + j3_z + j4_z + j5_z + j6_z
    ])
}

// ---------------------------- tests --------------------------------------
#[cfg(test)]
mod test {
    use super::*;
    // use all_asserts::assert_near;

    // #[test]
//     fn solar_position_test() {
//         let actual = arr1(&[1.43273853e+08, -4.64303244e+07, -2.10808621e+07]);

//         let function_output = solar_position(100.0);

//         assert_eq!(actual, function_output);
//     }

    #[test]
    fn test_exp_rho() {
        let val = exp_rho(500.0);

        let actual = 5.215e-13;

        assert_eq!(val, actual);
    }

    #[test]
    fn test_exp_rho_max() {
        let val = exp_rho(1000.0);

        let actual = 3.560997994537098e-15;

        assert_eq!(val, actual);
    }
}