use ark_ff::Field;
use std::ops::{Add, Mul, Neg, Sub};

/// Univariate polynomial with coefficients [c_0, c_1, ...] for c_0 + c_1x + ...
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnivariatePolynomial<F: Field> {
    pub coefficients: Vec<F>,
}

impl<F: Field> UnivariatePolynomial<F> {
    pub fn new(coefficients: &[F]) -> Self {
        Self {
            coefficients: coefficients.to_vec(),
        }
    }

    /// Evaluates the interpolated polynomial at a given point `x`.
    pub fn evaluate(&self, x: F) -> F {
        let mut result = F::zero();
        let mut x_power = F::one();

        for &coeff in &self.coefficients {
            result += coeff * x_power;
            x_power *= x;
        }

        result
    }

    /// Reduces the polynomial by removing trailing zeros.
    pub fn reduce(&mut self) {
        // Start from the last coefficient and remove zeros until a non-zero is found
        while let Some(&last) = self.coefficients.last() {
            if last == F::zero() {
                self.coefficients.pop();  // Remove the last element
            } else {
                break;  // Stop if the last element is non-zero
            }
        }
    }

    /// Checks if the polynomial is the zero polynomial (i.e., all coefficients are zero).
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(|&coeff| coeff == F::zero())
    }
}


// Implementing the Add trait for polynomial addition
impl<F: Field> Add for UnivariatePolynomial<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let max_len = self.coefficients.len().max(other.coefficients.len());
        let mut result_coeffs = vec![F::zero(); max_len];

        // Add coefficients from the first polynomial
        for i in 0..self.coefficients.len() {
            result_coeffs[i] += self.coefficients[i].clone();
        }

        // Add coefficients from the second polynomial
        for i in 0..other.coefficients.len() {
            result_coeffs[i] += other.coefficients[i].clone();
        }

        UnivariatePolynomial::new(&result_coeffs)
    }
}

// Implementing the Mul trait for polynomial multiplication
impl<F: Field> Mul for UnivariatePolynomial<F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut result_coeffs = vec![F::zero(); self.coefficients.len() + other.coefficients.len() - 1];

        // Perform the polynomial multiplication
        for i in 0..self.coefficients.len() {
            for j in 0..other.coefficients.len() {
                result_coeffs[i + j] += self.coefficients[i] * other.coefficients[j];
            }
        }

        UnivariatePolynomial::new(&result_coeffs)
    }
}

impl<F: Field> Neg for UnivariatePolynomial<F> {
    type Output = Self;

    fn neg(self) -> Self {
        // Negate each coefficient of the polynomial
        let neg_coeffs: Vec<F> = self.coefficients
            .into_iter()
            .map(|coeff| -coeff)
            .collect();

        UnivariatePolynomial::new(&neg_coeffs)
    }
}

impl<F: Field> Sub for UnivariatePolynomial<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let max_len = self.coefficients.len().max(other.coefficients.len());

        let mut result_coeffs = vec![F::zero(); max_len];

        // Subtract coefficients from the first polynomial
        for i in 0..self.coefficients.len() {
            result_coeffs[i] = result_coeffs[i] + self.coefficients[i].clone(); // Add first
        }

        // Subtract coefficients from the second polynomial
        for i in 0..other.coefficients.len() {
            result_coeffs[i] = result_coeffs[i] - other.coefficients[i].clone(); // Subtract second
        }

        UnivariatePolynomial::new(&result_coeffs)
    }
}


#[cfg(test)]
mod tests {
    use ark_ff::{Field, UniformRand, Zero};
    use ark_std::test_rng;
    use crate::constant_for_curves::ScalarField;
    use crate::ring::univariate_poly::UnivariatePolynomial;

    type F = ScalarField;

    /// Helper function to evaluate the polynomial P(x) = x^7 + 12 * x^6 + 5 * x^3 + 100 * x^2 + x + 9
    fn original_polynomial<F: Field>(x: F) -> F {
        x.pow([7]) +
            F::from(12u64) * x.pow([6]) +
            F::from(5u64) * x.pow([3]) +
            F::from(100u64) * x.pow([2]) +
            x +
            F::from(9u64)
    }

    #[test]
    fn univariate_evaluation_test() {
        let coeffs = [
            F::from(9u64), // x^0
            F::from(1u64), // x^1
            F::from(100u64), // x^2
            F::from(5u64), // x^3
            F::from(0u64), // x^4
            F::from(0u64), // x^5
            F::from(12u64), // x^6
            F::from(1u64), // x^7
        ];

        let poly = UnivariatePolynomial::new(&coeffs);
        let x = F::rand(&mut test_rng());

        // making sure the evaluation function works well
        assert_eq!(
            poly.evaluate(x),
            original_polynomial(x),
        );
    }

    #[test]
    fn univariate_multiplication_test() {
        // Represents 1 + 2x + 3x^2
        let p1 = UnivariatePolynomial::new(&[F::from(1u64), F::from(2u64), F::from(3u64)]);
        // Represents 4 + 5x
        let p2 = UnivariatePolynomial::new(&[F::from(4u64), F::from(5u64)]);

        assert_eq!(
            (p1 * p2).coefficients,
            vec![F::from(4u64), F::from(13u64), F::from(22u64), F::from(15u64)],
        );
    }

    #[test]
    fn univariate_reduction_test() {
        let mut p = UnivariatePolynomial::new(&[F::from(10u64), F::zero(), F::zero(), F::zero()]);
        p.reduce();

        assert_eq!(
            p.coefficients,
            vec![F::from(10u64)],
        );
    }

    #[test]
    fn test_sum_polynomials() {
        let p1 = UnivariatePolynomial::new(&[F::from(1u64), F::from(2u64), F::from(3u64)]);
        let p2 = UnivariatePolynomial::new(&[F::from(4u64), F::from(5u64)]);

        let sum_result = p1 + p2;
        assert_eq!(sum_result.coefficients, vec![F::from(5u64), F::from(7u64), F::from(3u64)]);
    }

    #[test]
    fn test_is_zero() {
        let p1 = UnivariatePolynomial::new(&[F::zero(), F::zero(), F::zero()]);
        let p2 = UnivariatePolynomial::new(&[F::from(1u64), F::zero(), F::from(2u64)]);

        assert!(p1.is_zero());  // Polynomial with all zero coefficients
        assert!(!p2.is_zero()); // Polynomial with non-zero coefficients
    }

    #[test]
    fn test_operator_behavior() {
        // Generate random coefficients for four polynomials
        let f1_coeffs: Vec<F> = (0..4).map(|_| F::rand(&mut test_rng())).collect();
        let f2_coeffs: Vec<F> = (0..4).map(|_| F::rand(&mut test_rng())).collect();
        let f3_coeffs: Vec<F> = (0..4).map(|_| F::rand(&mut test_rng())).collect();
        let f4_coeffs: Vec<F> = (0..4).map(|_| F::rand(&mut test_rng())).collect();

        // Create the polynomials
        let f1 = UnivariatePolynomial::new(&f1_coeffs);
        let f2 = UnivariatePolynomial::new(&f2_coeffs);
        let f3 = UnivariatePolynomial::new(&f3_coeffs);
        let f4 = UnivariatePolynomial::new(&f4_coeffs);

        // Random element r (here we'll use a random Fr value)
        let r = F::rand(&mut test_rng());

        // Evaluate the polynomials at r
        let f1_r = f1.evaluate(r);
        let f2_r = f2.evaluate(r);
        let f3_r = f3.evaluate(r);
        let f4_r = f4.evaluate(r);

        // Perform the operations (f1 + f2) * (f3 - f4) and evaluate at r
        let left_expr = (f1.clone() + f2.clone()) * (f3.clone() - f4.clone());
        let left_value = left_expr.evaluate(r);

        // Perform the operation [f1(r) + f2(r)] * [f3(r) - f4(r)]
        let right_value = (f1_r + f2_r) * (f3_r - f4_r);

        // Check that the results match
        assert_eq!(left_value, right_value);
    }
}