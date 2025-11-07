use ark_ff::PrimeField;
use ark_r1cs_std::{
    alloc::{AllocVar, AllocationMode},
    fields::fp::FpVar,
    prelude::*,
};
use ark_relations::r1cs::{ConstraintSystemRef, Namespace, SynthesisError};
use std::borrow::Borrow;
use ark_ec::pairing::Pairing;
use ark_std::{test_rng, UniformRand};
use crate::constant_for_curves::ScalarField;
use crate::merkle_tree::mimc::MiMC;

/// In-circuit MiMC variable version
pub struct MiMCVar<F: PrimeField> {
    pub rounds: usize,
    pub constants: Vec<FpVar<F>>,
}

impl<F: PrimeField> AllocVar<MiMC<F>, F> for MiMCVar<F> {
    fn new_variable<T: Borrow<MiMC<F>>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let binding = f()?;
        let mimc = binding.borrow();

        let constants = mimc
            .constants
            .iter()
            .map(|c| FpVar::new_variable(cs.clone(), || Ok(*c), mode))
            .collect::<Result<Vec<_>, _>>()?;

        Ok(MiMCVar {
            rounds: mimc.rounds,
            constants,
        })
    }
}

impl<F: PrimeField> MiMCVar<F> {
    /// In-circuit MiMC hash
    pub fn hash(&self, input: &FpVar<F>) -> Result<FpVar<F>, SynthesisError> {
        let mut x = input.clone();
        let k = FpVar::zero(); // fixed key = 0

        for i in 0..self.rounds {
            let c = &self.constants[i];
            let tmp = &x + &k + c;
            x = &tmp * &tmp * &tmp; // cube operation
        }

        Ok(x)
    }

    /// Sequential MiMC hashing (absorbing multiple inputs)
    pub fn hash_elements(&self, inputs: &[FpVar<F>]) -> Result<FpVar<F>, SynthesisError> {
        let mut acc = FpVar::zero();
        for x in inputs {
            acc += x;
            acc = self.hash(&acc)?;
        }
        Ok(acc)
    }
}

/// Generic MiMC test function over any pairing engine E
pub fn mr_path_verification<F: PrimeField, E: Pairing<ScalarField=F>>(
    cs: ConstraintSystemRef<F>,
    rounds: usize,
    n: usize,
) -> Result<(), SynthesisError> {
    // --- Step 1: Generate n random field elements ---
    let mut rng = test_rng();
    let inputs: Vec<F> = (0..n).map(|_| F::rand(&mut rng)).collect();

    let mimc = MiMC::<F>::new(rounds);
    let expected_hash = mimc.hash_elements(&inputs);

    // Allocate MiMC constants as constants
    let mimc_var = MiMCVar::new_variable(
        cs.clone(),
        || Ok(&mimc),
        AllocationMode::Constant,
    )?;

    // Allocate inputs as witnesses
    let input_vars = inputs
        .iter()
        .map(|x| FpVar::<F>::new_witness(cs.clone(), || Ok(*x)))
        .collect::<Result<Vec<_>, _>>()?;

    // --- Step 4: Compute MiMC hash inside the circuit ---
    let hash_var = mimc_var.hash_elements(&input_vars)?;

    // Allocate expected hash as a public input
    let expected_var = FpVar::<F>::new_input(cs.clone(), || Ok(expected_hash))?;

    // Enforce equality
    hash_var.enforce_equal(&expected_var)?;

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Fr;
    use ark_relations::r1cs::ConstraintSystem;

    #[test]
    fn test_mimc_var_consistency() {
        let rounds = 10;
        let inputs = vec![Fr::from(3u64), Fr::from(5u64), Fr::from(7u64)];

        // --- Outside-circuit computation ---
        let mimc = MiMC::<Fr>::new(rounds);
        let expected_hash = mimc.hash_elements(&inputs);

        // --- Build constraint system ---
        let cs = ConstraintSystem::<Fr>::new_ref();

        // Allocate MiMCVar in-constraint
        let mimc_var = MiMCVar::new_variable(
            cs.clone(),
            || Ok(&mimc),
            ark_r1cs_std::alloc::AllocationMode::Constant,
        )
            .unwrap();

        // Allocate inputs as witness variables
        let input_vars = inputs
            .iter()
            .map(|x| FpVar::new_witness(cs.clone(), || Ok(*x)))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        // Compute in-circuit hash
        let hash_var = mimc_var.hash_elements(&input_vars).unwrap();

        // Allocate expected output as public input
        let expected_var = FpVar::new_input(cs.clone(), || Ok(expected_hash)).unwrap();

        // Enforce equality
        hash_var.enforce_equal(&expected_var).unwrap();

        // --- Check constraint satisfaction ---
        assert!(cs.is_satisfied().unwrap());

        // --- Double-check concrete equality ---
        let evaluated_hash = hash_var.value().unwrap();
        assert_eq!(evaluated_hash, expected_hash);
    }
}
