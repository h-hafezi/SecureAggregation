use std::marker::PhantomData;
use ark_ff::{PrimeField};
use rand::RngCore;
use crate::ring::{Ring, RingParams};
use num::BigInt;

// Define the KAHE struct with parameters sigma_s and sigma_e
pub struct KAHE<ParamsQ: RingParams, ParamsT: RingParams> {
    _q: PhantomData<ParamsQ>,
    _t: PhantomData<ParamsT>,
    sigma_s: f64,  // Gaussian parameter for secret key distribution
    sigma_e: f64,  // Gaussian parameter for error distribution
}

impl<ParamsQ: RingParams, ParamsT: RingParams> KAHE<ParamsQ, ParamsT>
    where
        BigInt: From<<ParamsT as RingParams>::F>,
        BigInt: From<<ParamsQ as RingParams>::F> {
    // Constructor to initialize the KAHE struct with sigma_s and sigma_e
    pub fn new(sigma_s: f64, sigma_e: f64) -> Self {
        KAHE {
            _q: PhantomData,
            _t: PhantomData,
            sigma_s,
            sigma_e,
        }
    }

    // KAHE::Setup() - samples and returns a public parameter a ∈ R_{q}
    pub fn setup<R: ark_std::rand::RngCore>(&self, rng: &mut R) -> Ring<ParamsQ> {
        Ring::<ParamsQ>::rand(rng)
    }

    // KAHE.KeyGen() - samples and returns the secret key k ∈ R_{q} using sigma_s
    pub fn keygen<R: RngCore>(&self, rng: &mut R) -> Ring<ParamsQ> {
        Ring::<ParamsQ>::rand_gaussian(rng, self.sigma_s)
    }

    // KAHE.Enc(x, k) - samples error e from D_{σ_e} and encrypts the message x with key k
    pub fn encrypt<R: RngCore>(
        &self,
        rng: &mut R,
        a: Ring<ParamsQ>,
        x: Ring<ParamsT>,
        k: Ring<ParamsQ>,
    ) -> Ring<ParamsQ>
    where
        <<ParamsQ as RingParams>::F as PrimeField>::BigInt: From<BigInt>,
        <ParamsQ as RingParams>::F: From<<<ParamsT as RingParams>::F as PrimeField>::BigInt>
    {
        // Sample error e from D_{σ_e} (Gaussian with mean 0 and standard deviation sigma_e)
        let error = Ring::<ParamsQ>::rand_gaussian(rng, self.sigma_e);
        let t = <ParamsQ as RingParams>::F::from(<<ParamsT as RingParams>::F>::MODULUS);

        // c = a * k + t * e + x (mod q1)
        let mut c = a * k;  // a * k
        c = c + error.scalar_multiply(t);     // Add error e
        let x_q = Ring::<ParamsQ>::transform::<ParamsT>(&x);
        c = c + x_q;          // Add the message x

        c  // Return the ciphertext c
    }

    // KAHE.Dec(c, k) - decrypts the ciphertext c with the secret key k
    pub fn decrypt(
        &self,
        c: Ring<ParamsQ>,
        a: Ring<ParamsQ>,
        k: Ring<ParamsQ>,
    ) -> Ring<ParamsT>
    where
        <<ParamsT as RingParams>::F as PrimeField>::BigInt: From<BigInt>
    {
        // c - a * k (mod t1)
        let mut result = c.clone();
        result = result - (a * k);  // c - a * k
        let result = Ring::<ParamsT>::transform::<ParamsQ>(&result);

        result
    }
}
