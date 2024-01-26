use dotenv;
use ode_solvers::dop_shared::IntegrationError;
use std::fmt;
use std;
use std::num::ParseFloatError;
use std::str::ParseBoolError;

// create error handler for program
#[derive(Debug)]
pub enum PropagatorError {
    IntegratingError(String),
    FileWritingError(String),
    EnvVariableError(String),
    InputError(String),
}

// enable error to be displayed
impl fmt::Display for PropagatorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match self {
            PropagatorError::IntegratingError(string) => write!(f, "{}", string),
            PropagatorError::FileWritingError(string) => write!(f, "{}", string),
            PropagatorError::EnvVariableError(string) => write!(f, "{}", string),
            PropagatorError::InputError(string) => write!(f, "{}", string),
        }
    }
}

// transform various error types to internal error type
impl From<dotenv::Error> for PropagatorError {
    fn from(e: dotenv::Error) ->Self {
        PropagatorError::EnvVariableError(e.to_string())
    }
}

impl From<std::io::Error> for PropagatorError {
    fn from(e: std::io::Error) -> Self {
        PropagatorError::FileWritingError(e.to_string())
    }
}

impl From<IntegrationError> for PropagatorError {
    fn from(e: IntegrationError) -> Self {
        PropagatorError::IntegratingError(e.to_string())
    }
}

impl From<std::env::VarError> for PropagatorError {
    fn from(e: std::env::VarError) -> Self {
        PropagatorError::EnvVariableError(e.to_string())
    }
}

impl From<ParseFloatError> for PropagatorError {
    fn from(e: ParseFloatError) -> Self {
        PropagatorError::EnvVariableError(e.to_string())
    }
}

impl From<ParseBoolError> for PropagatorError {
    fn from(e: ParseBoolError) -> Self {
        PropagatorError::EnvVariableError(e.to_string())
    }
}