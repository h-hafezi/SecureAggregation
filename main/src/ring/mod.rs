use num::Zero;
use std::ops::{Add, Mul, Neg, Sub};
use ark_ff::Field;
use crate::univariate_poly::UnivariatePolynomial;

/// Define the Ring trait, which includes the field F and the value N.
pub trait RingParams {
    type F: Field;  // Associated field type
    const N: usize; // Associated constant N=2^k
}

/// the Ring structure
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring<Params: RingParams> {
    pub coefficients: Vec<Params::F>,
}

impl<Params: RingParams> Ring<Params> {
    pub fn new(coefficients: &[Params::F]) -> Self {
        assert_eq!(coefficients.len(), Params::N);
        Self { coefficients: coefficients.to_vec() }
    }

    // Scalar multiplication by a field element lambda.
    pub fn scalar_multiply(&self, lambda: Params::F) -> Self {
        let result: Vec<Params::F> = self.coefficients.iter()
            .map(|&coeff| coeff * lambda)
            .collect();

        Self::new(result.as_slice())
    }
}

impl<Params: RingParams> Add for Ring<Params> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let result = self.coefficients.iter()
            .zip(&other.coefficients)
            .map(|(a, b)| *a + *b)
            .collect::<Vec<Params::F>>();

        Self::new(result.as_slice())
    }
}

impl<Params: RingParams> Neg for Ring<Params> {
    type Output = Self;

    fn neg(self) -> Self {
        // Negate each coefficient of the polynomial
        let neg_coeffs: Vec<Params::F> = self.coefficients
                                                .into_iter()
                                                .map(|coeff| -coeff)
                                                .collect();

        Ring::new(neg_coeffs.as_slice())
    }
}

impl<Params: RingParams> Sub for Ring<Params> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl<Params: RingParams> Ring<Params> {
    pub fn reduce(univar: UnivariatePolynomial<Params::F>) -> Self {
        let length = univar.coefficients.len();

        if length < Params::N {
            // Case 1: If length is less than N, pad with zeros.
            Ring {
                coefficients: {
                    let mut padded_coeffs = univar.coefficients.clone();
                    padded_coeffs.resize(Params::N, Params::F::zero());
                    padded_coeffs
                },
            }
        } else if length == Params::N {
            // Case 2: If length is equal to N, return the polynomial as it is.
            Ring {
                coefficients: univar.coefficients,
            }
        } else if length <= 2 * Params::N {
            // Case 3: If length is greater than N but less than or equal to 2N, apply the reduction.
            let mut reduced_coeffs = vec![Params::F::zero(); Params::N];
            reduced_coeffs[..Params::N].copy_from_slice(&univar.coefficients[..Params::N]); // Efficient copy
            for i in Params::N..length {
                reduced_coeffs[i - Params::N] -= univar.coefficients[i].clone();
            }
            Self::new(reduced_coeffs.as_slice())
        } else {
            // Case for invalid length greater than or equal to 2N
            panic!("Polynomial length exceeds or equals 2N, which is not supported.")
        }
    }
}

#[cfg(test)]
mod test_reduction {
    use super::*;
    use ark_bn254::Fr as F;

    /// Example implementation of the Ring trait for a specific type.
    pub struct Params;
    impl RingParams for Params {
        type F = F;
        const N: usize = 8;
    }

    // Test case for when the length is less than N (padding is expected)
    #[test]
    fn test_reduce_short_polynomial() {
        let reduced = Ring::<Params>::reduce(
            UnivariatePolynomial {
            coefficients: vec![F::from(1), F::from(2)]
            }
        );

        // Expected output is a polynomial of length N with padding
        assert_eq!(reduced.coefficients.len(), Params::N);
        assert_eq!(reduced.coefficients[0], F::from(1));
        assert_eq!(reduced.coefficients[1], F::from(2));
        for i in 2..Params::N {
            assert_eq!(reduced.coefficients[i], F::zero());
        }
    }

    // Test case for when the length is exactly N (no modification)
    #[test]
    fn test_reduce_exact_polynomial() {
        let reduced = Ring::<Params>::reduce(
            UnivariatePolynomial {
            coefficients: (0..Params::N).map(|i| F::from(i as u64)).collect(),
            }
        );

        // No modification should happen, so the result should be the same as the input
        assert_eq!(reduced.coefficients.len(), Params::N);
        for i in 0..Params::N {
            assert_eq!(reduced.coefficients[i], F::from(i as u64));
        }
    }

    // Test case for when the length is greater than N but less than or equal to 2N (reduction occurs)
    #[test]
    fn test_reduce_long_polynomial() {
        let reduced = Ring::<Params>::reduce(
            UnivariatePolynomial {
                coefficients: (0..Params::N).map(|i| F::from(i as u64)).chain((0..Params::N).map(|i| F::from(i as u64))).collect(),
            }
        );

        // Expected result after reduction: first N elements should be the same,
        // and elements from N to 2N should be subtracted from their corresponding coefficients
        assert_eq!(reduced.coefficients.len(), Params::N);
        for i in 0..Params::N {
            assert_eq!(reduced.coefficients[i], F::zero());
        }
    }

    // Edge case: If the length is exactly 2N, ensure that we handle it properly
    #[test]
    #[should_panic(expected = "Polynomial length exceeds or equals 2N, which is not supported.")]
    fn test_reduce_too_large_polynomial() {
        let _ = Ring::<Params>::reduce(
            UnivariatePolynomial {
            coefficients: (0..2 * Params::N + 1).map(|i| F::from(i as u64)).collect(),
            }
        );
    }
}
