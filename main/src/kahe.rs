use std::marker::PhantomData;
use ark_ff::{BigInteger, PrimeField};
use rand::RngCore;
use crate::ring::{Ring, RingParams};

// Define the KAHE struct with parameters sigma_s and sigma_e
pub struct KAHE<ParamsQ: RingParams, ParamsT: RingParams> {
    _q: PhantomData<ParamsQ>,
    _t: PhantomData<ParamsT>,
    sigma_s: f64,  // Gaussian parameter for secret key distribution
    sigma_e: f64,  // Gaussian parameter for error distribution
}

impl<ParamsQ: RingParams, ParamsT: RingParams> KAHE<ParamsQ, ParamsT> {
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
    ) -> Ring<ParamsQ> {
        // Sample error e from D_{σ_e} (Gaussian with mean 0 and standard deviation sigma_e)
        let error = Ring::<ParamsQ>::rand_gaussian(rng, self.sigma_e);

        let t ={
            let modulus_bigint = <<ParamsT as RingParams>::F>::MODULUS;
            let modulus_bytes = modulus_bigint.to_bytes_le();

            // Now use the bytes to create an element of the desired field (ParamsQ::F)
            <ParamsQ as RingParams>::F::from_le_bytes_mod_order(modulus_bytes.as_slice())
        };

        // c = a * k + t * e + x (mod q1)
        let mut c = a * k;  // a * k
        c = c + error.scalar_multiply(t);     // Add error e
        let x_q = Ring::<ParamsQ>::transform::<ParamsT>(&x);
        c = c + x_q;          // Add the message x

        c
    }

    // KAHE.Dec(c, k) - decrypts the ciphertext c with the secret key k
    pub fn decrypt(
        &self,
        c: Ring<ParamsQ>,
        a: Ring<ParamsQ>,
        k: Ring<ParamsQ>,
    ) -> Ring<ParamsT> {
        // c - a * k (mod t1)
        let mut result = c.clone();
        result = result - (a * k);  // c - a * k
        let result = Ring::<ParamsT>::transform::<ParamsQ>(&result);

        result
    }
}

#[cfg(test)]
mod tests {
    use ark_bn254::Fr;
    use ark_ff::{Fp64, MontBackend, MontConfig};
    use ark_std::test_rng;
    use rand::{rng};
    use crate::kahe::KAHE;
    use crate::ring::{Ring, RingParams};

    // Define the configuration for Fp17 (a small prime field with p = 17)
    #[derive(MontConfig)]
    #[modulus = "17"]
    #[generator = "3"]
    pub struct Fp17Config;

    // Define Fp17 using a 64-bit Montgomery representation
    pub type Fp17 = Fp64<MontBackend<Fp17Config, 1>>;

    #[derive(Clone, Debug, PartialEq)]
    pub struct ParamsT;

    impl RingParams for ParamsT {
        type F = Fp17;
        const N: usize = 4;
    }

    // Define the large field ParamsQ using BN254 Fr
    #[derive(Clone, Debug)]
    pub struct ParamsQ;

    impl RingParams for ParamsQ {
        type F = Fr;
        const N: usize = 4;
    }

    #[test]
    fn test_kahe_encryption_decryption() {
        let sigma_s = 1.0;
        let sigma_e = 4.0;

        // Initialize the KAHE scheme
        let kahe = KAHE::<ParamsQ, ParamsT>::new(sigma_s, sigma_e);

        // Setup public parameters
        let a = kahe.setup(&mut test_rng());

        // Keygen to generate a secret key
        let k = kahe.keygen(&mut rng());

        // Encrypt a message (x) using the key (k)
        let x = Ring::<ParamsT>::rand(&mut test_rng());  // Sample message in ParamsT
        let ciphertext = kahe.encrypt(&mut rng(), a.clone(), x.clone(), k.clone());

        // Decrypt the message back
        let decrypted = kahe.decrypt(ciphertext, a, k);

        // Check if the decrypted message is the same as the original message
        assert_eq!(x, decrypted);
    }
}


