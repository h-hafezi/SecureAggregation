use num::{BigInt, Zero};
use std::ops::{Add, Mul, Neg, Sub};
use ark_ff::{Field, PrimeField, UniformRand};
use ark_serialize::CanonicalSerialize;
use rand::RngCore;
use crate::univariate_poly::UnivariatePolynomial;
use rand_distr::{Normal, Distribution};

/// Define the Ring trait, which includes the field F and the value N.
pub trait RingParams: Clone {
    type F: PrimeField;  // Associated field type
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

    /// Generates a random Ring element with coefficients of length Params::N
    pub fn rand<R: ark_std::rand::RngCore>(rng: &mut R) -> Self {
        let coefficients: Vec<Params::F> = (0..Params::N)
            .map(|_| Params::F::rand(rng)) // Generate random field elements
            .collect();

        Ring::<Params> { coefficients }
    }

    /// Generates a random Ring element with coefficients drawn from a Gaussian distribution
    /// with mean 0 and variance sigma^2 using a provided RNG.
    pub fn rand_gaussian<R: RngCore>(rng: &mut R, sigma: f64) -> Self {
        // Define a Gaussian distribution with mean 0 and variance sigma^2
        let normal_dist = Normal::new(0.0, sigma).unwrap();

        // Generate coefficients using the Gaussian distribution
        let coefficients: Vec<Params::F> = (0..Params::N)
            .map(|_| {
                // Sample from the Gaussian distribution
                let sample = normal_dist.sample(rng).round().abs() as u64;
                Params::F::from(sample)
            }).collect();

        Ring::<Params> { coefficients }
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

impl<ParamsF: RingParams> Ring<ParamsF> {
    pub fn transform<ParamsQ: RingParams>(ring: &Ring<ParamsQ>) -> Ring<ParamsF>
    where
        BigInt: From<<ParamsQ as RingParams>::F>,
        <<ParamsF as RingParams>::F as PrimeField>::BigInt: From<BigInt>,
    {
        let coefficients = ring.coefficients.iter().map(|x| {
            let value = BigInt::from(x.clone()); // Convert ParamsQ::F to BigInt

            // Convert BigInt to the correct F::BigInt representation
            let value_as_f_bigint = <<ParamsF as RingParams>::F as PrimeField>::BigInt::from(value);

            // Convert the BigInt representation into ParamsF::F
            ParamsF::F::from_bigint(value_as_f_bigint).unwrap()
        }).collect();

        Ring { coefficients }
    }
}

#[cfg(test)]
mod test {
    use ark_ff::{One, UniformRand, Zero};
    use ark_bn254::fr::Fr as F;
    use ark_std::test_rng;
    use crate::ring::{Ring, RingParams};

    #[derive(Clone)]
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

    #[test]
    fn test_ring_multiplication_identity() {
        let mut rng = test_rng();

        let f1 = Ring::<Params>::rand(&mut rng);
        let f2 = Ring::<Params>::rand(&mut rng);
        let f3 = Ring::<Params>::rand(&mut rng);
        let f4 = Ring::<Params>::rand(&mut rng);

        // Compute both sides of the identity:
        let lhs = (f1.clone() - f2.clone()) * (f3.clone() + f4.clone());
        let rhs = (f1.clone() * f3.clone()) - (f2.clone() * f3.clone()) - (f2.clone() * f4.clone()) + (f1.clone() * f4.clone());

        // Check if both sides are equal
        assert_eq!(
            lhs.coefficients,
            rhs.coefficients,
            "Multiplication identity (f1 - f2) * (f3 + f4) != (f1 * f3) - (f2 * f3) - (f2 * f4) + (f1 * f4)"
        );
    }
}

