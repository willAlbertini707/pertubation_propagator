use ode_solvers::*;

use dotenv::dotenv;
use orbit_tools::errors::PropagatorError;
use std::env;
use std::path::Path;

// internal imports
mod orbit_tools;
use orbit_tools::propagator::OrbitPropagator;

type State = Vector6<f64>;
fn main() {
    
    // try to run propagator 
    match run_propagator_from_env() {
        Ok(_) => println!("Propagation was successfully run"),
        Err(e) => println!("An error was encountered: {:?}", e)
    }
}

fn run_propagator_from_env() -> Result<(), PropagatorError>{

    // load env variables into scope
    dotenv()?;

    // pull out states
    let initial_state_string = env::var("initial_state")?;
    let initial_state = string_list_to_vector(initial_state_string)?;
    let area = env::var("area")?.parse::<f64>()?;
    let cd = env::var("cd")?.parse::<f64>()?;
    let mass = env::var("mass")?.parse::<f64>()?;
    let mu = env::var("mu")?.parse::<f64>()?;
    let t0 = env::var("t0")?.parse::<f64>()?;
    let tf = env::var("tf")?.parse::<f64>()?;
    let dx = env::var("dx")?.parse::<f64>()?;
    let init_jd = env::var("init_jd")?.parse::<f64>()?;
    let cr = env::var("cr")?.parse::<f64>()?;
    let body_radius = env::var("body_radius")?.parse::<f64>()?;
    let add_pertibutaions = env::var("add_pertubations")?.parse::<bool>()?;

    // create var with long enough lifetime
    let file_name = env::var("output_file")?;
    let output_file = Path::new(&file_name[..]);


    // propogate orbit based on values
    let mut my_orbit = OrbitPropagator::new(initial_state, 
                                                                area, 
                                                                cd, 
                                                                mass, 
                                                                mu, 
                                                                init_jd, 
                                                                cr, 
                                                                body_radius,
                                                                add_pertibutaions
                                                            );
    my_orbit.propogate(t0, tf, dx)?;
    
    my_orbit.write_to_file(output_file)?;
    Ok(())
}

fn string_list_to_vector(string_list: String) -> Result<State, PropagatorError> {
    // helper function to read in initial state and convert
    // from string to vector

    let values: Vec<f64> = string_list.split(",").map(|a| a.parse::<f64>().unwrap()).collect();

    // make sure initial state only has 6 values
    if values.len() != 6 {
        return Err(PropagatorError::InputError("State data not right size".to_string()))
    }

    Ok(Vector6::from_vec(values))
}