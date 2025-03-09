use crate::univariate_poly::UnivariatePolynomial;
use ark_ff::{BigInteger, PrimeField, UniformRand};
use num::Zero;
use rand::{Rng, RngCore};
use rand_distr::{Distribution, Normal};
use std::ops::{Add, Mul, Neg, Sub};

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
    pub fn rand_gaussian<R: RngCore + Rng>(rng: &mut R, sigma: f64) -> Self {
        // Define a Gaussian distribution with mean 0 and variance sigma^2
        let normal_dist = Normal::new(0.0, sigma).unwrap();

        // Generate coefficients using the Gaussian distribution
        let coefficients: Vec<Params::F> = (0..Params::N)
            .map(|_| {
                // Sample from the Gaussian distribution
                let sample = normal_dist.sample(rng) as u64;
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

impl<ParamsF: RingParams> Ring<ParamsF>  {
    pub fn transform<ParamsQ: RingParams>(ring: &Ring<ParamsQ>) -> Ring<ParamsF> {
        assert_eq!(ParamsF::N, ParamsQ::N, "cannot do conversion between these rings");

        let coefficients = ring.coefficients.iter().map(|x| {
            // Use the cast_field function to convert from ParamsQ::F to ParamsF::F
            cast_field::<ParamsQ::F, ParamsF::F>(*x)
        }).collect();

        Ring { coefficients }
    }
}

pub fn cast_field<Fr, Fq>(first_field: Fr) -> Fq
where
    Fr: PrimeField,
    Fq: PrimeField,
{
    // Convert the Fr element to its big integer representation
    let bytes = first_field.into_bigint().to_bytes_le();

    // Convert the big integer representation to an Fq element
    let fq_element = Fq::from_le_bytes_mod_order(bytes.as_slice());

    fq_element
}

#[cfg(test)]
mod test {
    use crate::constant_curve::ScalarField as F;
    use crate::ring::{Ring, RingParams};
    use ark_ff::{One, Zero};
    use ark_std::test_rng;

    #[derive(Clone)]
    pub struct Params;

    impl RingParams for Params {
        type F = F;
        const N: usize = 4;
    }

    #[test]
    fn test_ring_multiplication() {
        let one = F::one();
        let two = &one + &one;

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

        let expected = Ring::<Params> { coefficients: vec![two, F::zero(), two, F::zero()] };

        let result = f * g;

        assert!(
            result.coefficients.iter().zip(expected.coefficients.iter()).all(|(a, b)| *a == *b),
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
        assert!(
            lhs.coefficients.iter().zip(rhs.coefficients.iter()).all(|(a, b)| *a == *b),
            "Multiplication identity (f1 - f2) * (f3 + f4) != (f1 * f3) - (f2 * f3) - (f2 * f4) + (f1 * f4)"
        );
    }
}

