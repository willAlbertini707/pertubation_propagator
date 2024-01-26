// The equations of motion describing the motion of a spacecraft on a Kepler
// orbit are integrated using Rk4.
use ode_solvers::rk4::*;
use ode_solvers::*;
use std::path::Path;
use std::{fs::File, io::BufWriter, io::Write};

use ndarray::{Array1, arr1};

// internal exports
use super::errors::PropagatorError;
use super::pertubations::{drag_acceleration, three_body_solar, solar_position, srp, j2_j6};

type State = Vector6<f64>;
type Time = f64;

// system object
#[derive(Clone)]
pub struct OrbitSystem {
    pub mu: f64,
    pub area: f64,
    pub cd: f64,
    pub mass: f64,
    pub init_jd: f64,
    pub cr: f64,
    pub body_radius: f64,
    pub add_petubations: bool,
}

// propagator object
pub struct OrbitPropagator {
    init_state: State,
    state_data: Option<Vec<State>>,
    time: Option<Vec<Time>>,
    body: OrbitSystem,
}

impl OrbitPropagator {

    pub fn new(init_state: State, 
                area: f64, 
                cd: f64, 
                mass: f64,
                mu: f64, 
                init_jd: f64, 
                cr: f64, 
                body_radius: f64,
                add_pertubations: bool
            ) -> Self {
        // create the OrbitPropagator object

        OrbitPropagator { 
            init_state: init_state, 
            state_data: None, 
            time: None, 
            body: OrbitSystem{
                mu: mu,
                area: area,
                cd: cd,
                mass: mass,
                init_jd: init_jd,
                cr: cr,
                body_radius: body_radius,
                add_petubations: add_pertubations,

            } }
    }

    pub fn propogate(&mut self, t0: f64, tf: f64, dx: f64) -> Result<(), PropagatorError> {
        // propogate the orbit and store results in self

        // use Rk4 for propogation
        let mut stepper = Rk4::new(
            self.body.clone(),
            t0,
            self.init_state,
            tf,
            dx
        );
        let _res = stepper.integrate()?;

        // return independent variable
        self.time = Some(stepper.x_out().to_owned());

        // return state
        self.state_data = Some(stepper.y_out().to_owned());

        Ok(())
    }


    pub fn write_to_file(&self, file_path: &Path) -> Result<(), PropagatorError> {
        // write the results to a txt file

        if let None = self.state_data {
            // make sure orbit is propogated first
            println!("Propogate must be run before file can be written")
        }
        
        // create file
        let file = File::create(file_path)?;

        // create buf writer
        let mut buf = BufWriter::new(file);

        // write data in csv format
        for (i, state) in self.state_data.clone().unwrap().iter().enumerate() {

            buf.write_fmt(format_args!("{}", self.time.as_ref().unwrap()[i]))?;

            for val in state.iter() {
                buf.write_fmt(format_args!(", {}", val))?;
            }
            buf.write_fmt(format_args!("\n"))?;
        }

        buf.flush()?;

        Ok(())
    }
}

// implement System/equations for system data object
impl ode_solvers::System<State> for OrbitSystem {

    fn system(&self, t:Time, y: &State, dy: &mut State) {

        let r= arr1(&y.as_slice()[..3]);
        let v = arr1(&y.as_slice()[3..]);

        let r_mag = r.dot(&r).sqrt();

        let a_drag = drag_acceleration(r.clone(), r_mag, v.clone(), &self);

        // calculate current sun position
        let current_time = self.init_jd + t / 3600.0 / 24.0;
        let r_sun = solar_position(current_time);
        let r_sun_mag = r_sun.dot(&r_sun).sqrt();

        // calculate solar n-body
        let n_accel = three_body_solar(r.clone(), r_sun.clone(), r_sun_mag);

        // calculate SRP
        let srp_accel = srp(r.clone(), r_mag, r_sun.clone(), r_sun_mag, &self);

        // calculate J2-J6
        let j_accel = j2_j6(r.clone(), r_mag, &self);

        dy[0] = y[3];
        dy[1] = y[4];
        dy[2] = y[5];

        // add accelerations

        let mut ddy = -self.mu * r / r_mag.powi(3);

        if self.add_petubations {
            ddy = ddy +  a_drag + n_accel + srp_accel + j_accel;
        }

        // unpack into state vector
        dy[3] = ddy[0];
        dy[4] = ddy[1];
        dy[5] = ddy[2];

    }

    fn solout(&mut self, _t: Time, y: &State, _dy: &State) -> bool {
        true
    }
}