use core::mem;

/// A naive implementation of the Runge-Kutta explicit Runge-Kutta method.
///
/// Each step is computed following the formula:
///
/// ```math
/// y_{n+1} = y_n + h \sum_{i=1}^s b_i k_i
/// ```
///
/// where `k_i` is the `i`-th stage of the method, computed with
///
/// ```math
/// k_i = f(t_n + c_i h, y_n + h \sum_{j=1}^{i-1} a_{ij} k_j)
/// ```
#[derive(Debug)]
pub struct Naive<T, const S: usize> {
    /// The coefficients `$a_{ij}$` of the Runge-Kutta method.
    pub matrix: [[T; S]; S],
    /// The vector of weights `$b_i$` of the Runge-Kutta method.
    pub weights: [T; S],
    /// The vector of nodes `$c_i$` of the Runge-Kutta method.
    pub nodes: [T; S],
}

impl<T, const S: usize> Naive<T, S>
where
    T: num::Zero,
{
    /// Creates a new instance of the method.
    ///
    /// # Errors
    ///
    /// This will perform some checks on the values to ensure that they are
    /// consistent with a naive implementation of the explicit Runge-Kutta
    /// method. In particular, the following checks will be performed:
    ///
    /// - The matrix `$a_{ij}$` must be lower triangular, with the diagonal
    ///   elements equal to zero.
    ///
    /// - The matrix `$a_{ij}$` must be at least of size `$s \cross s$`.
    ///
    /// - The vector of weights and nodes must be of length `$s$`.
    ///
    /// Note that in all of the above.
    pub fn new(
        matrix: impl IntoIterator<Item = impl IntoIterator<Item = T>>,
        weights: impl IntoIterator<Item = T>,
        nodes: impl IntoIterator<Item = T>,
    ) -> Result<Self, NaiveError> {
        let weights = <[T; S]>::try_from(weights.into_iter().collect::<Vec<T>>())
            .map_err(|_| NaiveError::WeightsDim)?;
        let nodes = <[T; S]>::try_from(nodes.into_iter().collect::<Vec<T>>())
            .map_err(|_| NaiveError::NodesDim)?;

        // Unlike the other two attributes, the matrix is a bit more involved in
        // initialising due to Rust's strict type system. We have to construct
        // the slice which leaves it temporarily uninitialised.
        let matrix = {
            let mut tmp = mem::MaybeUninit::<[[T; S]; S]>::uninit();
            let ptr = tmp.as_mut_ptr() as *mut [[T; S]; S];

            let mut rows = matrix.into_iter();
            for i in 0..S {
                let row: Vec<_> = rows
                    .next()
                    .ok_or(NaiveError::MatrixDim)?
                    .into_iter()
                    .collect();
                let tmp_row = <[T; S]>::try_from(row).map_err(|_| NaiveError::MatrixDim)?;
                unsafe {
                    (*ptr)[i] = tmp_row;
                }
            }

            unsafe { tmp.assume_init() }
        };

        let naive = Self {
            matrix,
            weights,
            nodes,
        };

        // Ensure that the matrix is strictly lower triangular.
        for (i, row) in naive.matrix.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                if i < j && !value.is_zero() {
                    return Err(NaiveError::NonLowerTriangularMatrix);
                }
            }
        }

        Ok(naive)
    }
}

#[derive(Debug)]
pub enum NaiveError {
    MatrixDim,
    WeightsDim,
    NodesDim,
    NonLowerTriangularMatrix,
}

impl std::fmt::Display for NaiveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MatrixDim => {
                write!(f, "The matrix has the wrong dimension.")
            }
            Self::WeightsDim => {
                write!(f, "The weights vector has the wrong dimension.")
            }
            Self::NodesDim => {
                write!(f, "The nodes vector has the wrong dimension.")
            }
            Self::NonLowerTriangularMatrix => {
                write!(f, "The matrix is not strictly lower triangular.")
            }
        }
    }
}

impl std::error::Error for NaiveError {}

#[cfg(test)]
mod tests {
    use std::error;

    #[test]
    fn naive_new() -> Result<(), Box<dyn error::Error>> {
        super::Naive::<_, 2>::new(
            vec![vec![0.0, 0.0], vec![0.5, 0.0]],
            vec![0.0, 1.0],
            vec![0.0, 0.5],
        )?;

        Ok(())
    }

    #[test]
    fn naive_new_non_triangular() -> Result<(), Box<dyn error::Error>> {
        match super::Naive::<_, 2>::new(
            vec![vec![1.0, 2.0], vec![3.0, 4.0]],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        ) {
            Err(super::NaiveError::NonLowerTriangularMatrix) => (),
            Err(e) => Err(format!("Expected `NonLowerTriangularMatrix`, got {:?}", e))?,
            Ok(_) => Err("Expected `NonLowerTriangularMatrix`, got `Ok`")?,
        }

        Ok(())
    }

    #[test]
    fn naive_new_matrix_wrong_dim() -> Result<(), Box<dyn error::Error>> {
        match super::Naive::<_, 2>::new(
            vec![vec![0.0], vec![1.0, 0.0]],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        ) {
            Err(super::NaiveError::MatrixDim) => (),
            Err(e) => Err(format!("Expected `MatrixDim`, got {:?}", e))?,
            Ok(_) => Err("Expected `MatrixDim`, got `Ok`")?,
        };

        match super::Naive::<_, 2>::new(
            vec![vec![], vec![1.0, 0.0]],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        ) {
            Err(super::NaiveError::MatrixDim) => (),
            Err(e) => Err(format!("Expected `MatrixDim`, got {:?}", e))?,
            Ok(_) => Err("Expected `MatrixDim`, got `Ok`")?,
        }

        Ok(())
    }

    #[test]
    fn naive_new_weights_wrong_dim() -> Result<(), Box<dyn error::Error>> {
        match super::Naive::<_, 2>::new(
            vec![vec![0.0, 0.0], vec![1.0, 0.0]],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0],
        ) {
            Err(super::NaiveError::WeightsDim) => Ok(()),
            Err(e) => Err(format!("Expected `WeightsDim`, got {:?}", e))?,
            Ok(_) => Err("Expected `WeightsDim`, got `Ok`")?,
        }

        match super::Naive::<_, 2>::new(
            vec![vec![0.0, 0.0], vec![1.0, 0.0]],
            vec![0.0],
            vec![0.0, 0.0],
        ) {
            Err(super::NaiveError::WeightsDim) => Ok(()),
            Err(e) => Err(format!("Expected `WeightsDim`, got {:?}", e))?,
            Ok(_) => Err("Expected `WeightsDim`, got `Ok`")?,
        }
    }

    #[test]
    fn naive_new_nodes_wrong_dim() -> Result<(), Box<dyn error::Error>> {
        match super::Naive::<_, 2>::new(
            vec![vec![0.0, 0.0], vec![1.0, 0.0]],
            vec![0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ) {
            Err(super::NaiveError::NodesDim) => (),
        }
        match super::Naive::<_, 2>::new(
            vec![vec![0.0, 0.0], vec![1.0, 0.0]],
            vec![0.0, 0.0],
            vec![0.0],
        ) {
            Err(super::NaiveError::NodesDim) => (),
        }
    }
}
