//! Define the traits to handle initial value problems and the corresponding
//! solvers.
use core::fmt;
use std::error;

use crate::system::System;

/// Initial value problems (IVP), where the value of the function is known at a
/// point in time, and the goal is to find the value of the function at other
/// times. Despite it's name, the known point need not be the initial value.
///
/// The initial value problem is defined by
///
/// ```math
/// y(t_0) = y_0
/// ```
pub trait InitialValueProblem<T, Y>: System<T, Y> {
    /// Specify the initial value for the system to be solved.
    fn initial_value(&mut self, t0: T, y0: Y);
}

#[derive(Debug)]
pub enum Error {
    /// The solver failed to converge.
    ConvergenceFailed,
    /// The solver failed to converge within the specified number of iterations.
    MaxIterationsExceeded,
    /// The solver failed to converge within the specified tolerance.
    ToleranceExceeded,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::ConvergenceFailed => write!(f, "Convergence failed"),
            Error::MaxIterationsExceeded => write!(f, "Maximum number of iterations exceeded"),
            Error::ToleranceExceeded => write!(f, "Tolerance exceeded"),
        }
    }
}

impl error::Error for Error {}

pub trait SolverBuilder<T, Y> {
    type Solver: Solver<T, Y>;

    /// Build the solver.
    fn build(self) -> Self::Solver;
}

/// Generic solver for an initial value problem.
///
/// Implementation of this will typically be instantiated through a
/// corresponding `SolverBuilder`, and each solver will typically have a unique
/// set of options to control the algorithm.
pub trait Solver<T, Y> {
    /// Step the system by `dt`.
    ///
    /// The value of `dt` may be negative to step backwards in time.
    fn step(&mut self, dt: T) -> Y;

    /// Solve the system until time `$t$` is reached.
    ///
    /// The value of `t` may be less than the current time in which case the
    /// solver will step backwards until the time is reached.
    ///
    /// The parameters used internally when solving the system
    fn solve(&mut self, t: T) -> Result<Y, Error>;
}

/// Embedded solvers for initial value problems.
pub trait EmbeddedSolver<T, Y> {
    /// Get the error estimate for the last step performed.
    ///
    /// The error estimate is the difference between the solution obtained with
    /// the embedded solver and the solution obtained with the main solver for
    /// the last step. It is given by:
    ///
    /// ```math
    /// e_i \defeq h \sum_{i=j}^s d_j k_j
    /// ```
    ///
    /// where `$h$` is the last step size and `$d_j \defeq b_j^* - b_j$` is the
    /// difference between the weights of the main method and the embedded
    /// method.
    fn error_estimate(&self) -> Y;

    /// Compute the next step size based on the error estimate.
    fn step_size(&self) -> T;
}
