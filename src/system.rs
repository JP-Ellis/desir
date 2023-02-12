//! Generic traits to define differential equations.
//!
//! Unless otherwise specified, it is assumed that the differentiatial equation
//! is specified by the equation
//!
//! ```math
//! \ddfrac{y}{t} = f(t, y)
//! ```
//!
//! where `$y$` will, in general, by a vector-valued function.
//!
//! The trait [`System`] is used to define the differential equation which is to
//! be solved, and the library defines a number of solvers for the various
//! problems one might wish to solve relating to a particular system.

/// Generic trait for a differential equation system.
///
/// The differential equation system is defined by
///
/// ```math
/// \ddfrac{y}{t} = f(t, y)
/// ```
///
/// The time variable `$t$` is of type `T` and the state `$y$` is of type `Y`.
pub trait System<T, Y> {
    /// Compute the derivative of the state `$y$` with respect to `$t$` at the
    /// specified time.
    ///
    /// # Implementation
    ///
    /// The function is designed to consume the state vector `$y$` and return a
    /// value of the same type. This allows the implementation to reuse the
    /// memory allocated for the state vector, which can be useful for
    /// performance reasons.
    fn eval(&mut self, t: &T, y: Y) -> Y;
}
