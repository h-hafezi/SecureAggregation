use ark_ff::PrimeField;

/// Simple MiMC hash for field elements
pub struct MiMC<F: PrimeField> {
    pub(crate) rounds: usize,
    pub(crate) constants: Vec<F>,
}

impl<F: PrimeField> MiMC<F> {
    /// Create new MiMC instance with given number of rounds
    pub fn new(rounds: usize) -> Self {
        // Use simple constants: 1,2,3,... or can be from a hash
        let constants = (0..rounds)
            .map(|i| F::from(i as u64))
            .collect::<Vec<_>>();
        Self { rounds, constants }
    }

    /// Hash a single field element (can be used for leaves)
    pub fn hash(&self, input: F) -> F {
        let mut x = input;
        let k = F::zero(); // key = 0
        for i in 0..self.rounds {
            let c = self.constants[i];
            x = (x + k + c).pow([3u64]); // MiMC exponent = 3
        }
        x
    }

    /// Hash multiple field elements sequentially (absorbing)
    pub fn hash_elements(&self, inputs: &[F]) -> F {
        let mut acc = F::zero();
        for x in inputs {
            acc += *x;
            acc = self.hash(acc);
        }
        acc
    }
}
