use num::Zero;
use std::ops::{Add, Mul};
use ark_ff::Field;
use crate::univariate_poly::UnivariatePolynomial;

const N: usize = 8;  // Example N = 8, you can adjust as needed

// Define the RingPolynomial structure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingPolynomial<F: Field> {
    pub coefficients: Vec<F>,  // Coefficients of the polynomial (c_0, c_1, ..., c_{N-1}).
}

impl<F: Field> RingPolynomial<F> {
    // Constructor for a new RingPolynomial with coefficients.
    pub fn new(coefficients: Vec<F>) -> Self {
        assert_eq!(coefficients.len(), N);  // Ensure that the size is fixed at N.
        Self { coefficients }
    }

    // Polynomial addition (coefficient-wise).
    pub fn add(&self, other: &Self) -> Self {
        let mut result = Vec::with_capacity(self.coefficients.len());
        for (a, b) in self.coefficients.iter().zip(&other.coefficients) {
            result.push(*a + *b);
        }
        Self::new(result)
    }

    // Scalar multiplication by a field element lambda.
    pub fn scalar_multiply(&self, lambda: F) -> Self {
        let result: Vec<F> = self.coefficients.iter()
            .map(|&coeff| coeff * lambda.clone())
            .collect();
        Self::new(result)
    }

    // Negacyclic multiplication using matrix-vector multiplication.
    pub fn multiply(&self, other: &Self) -> Self {
        let mut result = vec![F::zero(); N];
        // Construct the negacyclic matrix for `self`.
        let negacyclic_self = self.negacyclic_matrix();

        // Perform the matrix-vector multiplication.
        for i in 0..N {
            for j in 0..N {
                result[i] = result[i].clone() + negacyclic_self[i][j].clone() * other.coefficients[j].clone();
            }
        }

        // Return the product polynomial
        Self::new(result)
    }

    // Helper function to construct the negacyclic matrix for the polynomial.
    fn negacyclic_matrix(&self) -> Vec<Vec<F>> {
        let mut matrix = vec![vec![F::zero(); N]; N];
        for i in 0..N {
            for j in 0..N {
                let idx = (i + j) % N;
                matrix[i][j] = self.coefficients[idx].clone();
            }
        }
        matrix
    }
}

impl<F: Field> RingPolynomial<F> {
    pub fn reduce(univar: UnivariatePolynomial<F>) -> Self {
        let length = univar.coefficients.len();
        match length {
            // Case 1: If length is less than N, pad with zeros.
            len if len < N => RingPolynomial {
                coefficients: {
                    let mut padded_coeffs = univar.coefficients.clone();
                    padded_coeffs.resize(N, F::zero());
                    padded_coeffs
                },
            },
            // Case 2: If length is equal to N, return the polynomial as it is.
            N => RingPolynomial {
                coefficients: univar.coefficients,
            },
            // Case 3: If length is greater than N but less than or equal to 2N, apply the reduction.
            len if len <= 2 * N => {
                let mut reduced_coeffs = vec![F::zero(); N];
                reduced_coeffs[..N].copy_from_slice(&univar.coefficients[..N]); // Efficient copy
                for i in N..len {
                    reduced_coeffs[i - N] -= univar.coefficients[i].clone();
                }
                Self::new(reduced_coeffs)
            },
            // Case for invalid length greater than or equal to 2N
            _ => panic!("Polynomial length exceeds or equals 2N, which is not supported."),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr as F;


    // Test case for when the length is less than N (padding is expected)
    #[test]
    fn test_reduce_short_polynomial() {
        let univar = UnivariatePolynomial {
            coefficients: vec![F::from(1), F::from(2)]
        };
        let reduced = RingPolynomial::<F>::reduce(univar);

        // Expected output is a polynomial of length N with padding
        assert_eq!(reduced.coefficients.len(), N);
        assert_eq!(reduced.coefficients[0], F::from(1));
        assert_eq!(reduced.coefficients[1], F::from(2));
        for i in 2..N {
            assert_eq!(reduced.coefficients[i], F::zero());
        }
    }

    // Test case for when the length is exactly N (no modification)
    #[test]
    fn test_reduce_exact_polynomial() {
        let univar = UnivariatePolynomial {
            coefficients: (0..N).map(|i| F::from(i as u64)).collect(),
        };
        let reduced = RingPolynomial::<F>::reduce(univar);

        // No modification should happen, so the result should be the same as the input
        assert_eq!(reduced.coefficients.len(), N);
        for i in 0..N {
            assert_eq!(reduced.coefficients[i], F::from(i as u64));
        }
    }

    // Test case for when the length is greater than N but less than or equal to 2N (reduction occurs)
    #[test]
    fn test_reduce_long_polynomial() {
        let univar = UnivariatePolynomial {
            coefficients: (0..N).map(|i| F::from(i as u64)).chain((0..N).map(|i| F::from(i as u64))).collect(),
        };
        let reduced = RingPolynomial::<F>::reduce(univar);

        // Expected result after reduction: first N elements should be the same,
        // and elements from N to 2N should be subtracted from their corresponding coefficients
        assert_eq!(reduced.coefficients.len(), N);
        for i in 0..N {
            assert_eq!(reduced.coefficients[i], F::zero());
        }
    }

    // Edge case: If the length is exactly 2N, ensure that we handle it properly
    #[test]
    #[should_panic(expected = "Polynomial length exceeds or equals 2N, which is not supported.")]
    fn test_reduce_too_large_polynomial() {
        let univar = UnivariatePolynomial {
            coefficients: (0..2 * N + 1).map(|i| F::from(i as u64)).collect(),
        };
        let _ = RingPolynomial::<F>::reduce(univar);
    }
}
