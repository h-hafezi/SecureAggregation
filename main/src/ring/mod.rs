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

    /// Checks if all coefficients are zero.
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|&c| c == Params::F::zero())
    }

    /// Scalar multiplication by a field element lambda.
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

impl<Params: RingParams> Mul for Ring<Params> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let poly1 = UnivariatePolynomial {
            coefficients: self.coefficients.clone(),
        };

        let poly2 = UnivariatePolynomial {
            coefficients: other.coefficients.clone(),
        };

        let poly_mult = poly1 * poly2;
        assert_eq!(poly_mult.coefficients.len(), 2 * Params::N - 1);

        let mut reduced_coeffs = vec![Params::F::zero(); Params::N];
        reduced_coeffs[..Params::N].copy_from_slice(&poly_mult.coefficients[..Params::N]); // Efficient copy
        for i in Params::N..2 * Params::N - 1 {
            reduced_coeffs[i - Params::N] -= poly_mult.coefficients[i].clone();
        }

        Ring {
            coefficients: reduced_coeffs,
        }
    }
}

#[cfg(test)]
mod test {
    use ark_ff::{One, Zero};
    use ark_bn254::fr::Fr as F;
    use crate::ring::{Ring, RingParams};

    pub struct Params;

    impl RingParams for Params {
        type F = F;
        const N: usize = 4;
    }

    #[test]
    fn test_ring_multiplication() {
        let one = F::one();

        // Define two polynomials of degree 3:
        // f(x) = 1 + x + x^2 + x^3
        // g(x) = 1 - x + x^2 - x^3
        let f = Ring::<Params> { coefficients: vec![one, one, one, one] };
        let g = Ring::<Params> { coefficients: vec![one, -one, one, -one] };

        // Expected result before reduction:
        // f(x) * g(x) = (1 + x + x^2 + x^3) * (1 - x + x^2 - x^3)
        //             = (1 - x + x^2 - x^3) + (x - x^2 + x^3 - x^4)
        //             + (x^2 - x^3 + x^4 - x^5) + (x^3 - x^4 + x^5 - x^6)
        // => 1 + 0x + 1x^2 + 0x^3 - 1x^4 - 0x^5 - 1x^6
        // => Reduced result: 2 + 0x + 2x^2 + 0x^3

        let expected = Ring::<Params> { coefficients: vec![F::from(2), F::zero(), F::from(2), F::zero()] };

        let result = f * g;

        assert_eq!(
            result.coefficients,
            expected.coefficients,
            "Multiplication result is incorrect."
        );
    }

    #[test]
    fn test_ring_multiplication_different_polynomials() {
        let one = F::one();
        let two = one + one;
        let three = two + one;
        let four = two + two;

        // Define two polynomials of degree 3:
        // h(x) = 1 + 2x + 3x^2 + 4x^3
        // k(x) = 4 + 3x + 2x^2 + x^3
        let h = Ring::<Params> { coefficients: vec![one, two, three, four] };
        let k = Ring::<Params> { coefficients: vec![four, three, two, one] };

        // Expected result after reduction:
        let expected = Ring::<Params> { coefficients: vec![-four * four, F::zero(), four * four, F::from(30)] };

        let result = h * k;

        assert_eq!(
            result.coefficients,
            expected.coefficients,
            "Multiplication result is incorrect."
        );
    }

}

