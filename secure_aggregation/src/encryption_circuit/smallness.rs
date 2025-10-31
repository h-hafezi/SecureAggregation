use ark_ff::{BigInteger, PrimeField};
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::prelude::*;
use ark_relations::r1cs::SynthesisError;
use ark_std::vec::Vec;
use std::ops::Mul;
use crate::commitment::Len;

/// Enforce that `x_var < limit` by decomposing `x_var` into its binary bits.
///
/// Since `limit` is assumed to be a power of two (`limit = 2^k`), we can
/// represent `x_var` as a sum of its bits multiplied by powers of two:
///
///     x_var = sum_{i=0}^{k-1} 2^i * bit_i
///
/// Each `bit_i` is constrained to be boolean, which guarantees
/// that 0 <= x_var < 2^k.
///
/// # Arguments
/// - `x_var`: the value to constrain
/// - `limit`: must equal 2^k (constant in the circuit)
/// - `pow2_coeffs`: precomputed powers of two [1, 2, 4, ..., 2^(k-1)] as FpVar constants
pub fn enforce_x_less_than_limit_power_of_two<F: PrimeField>(
    x_var: &FpVar<F>,
    n_bits: usize,
    pow2_coeffs: &[FpVar<F>],
) -> Result<(), SynthesisError> {
    // assert the length of pow2_coeffs is correct
    assert_eq!(pow2_coeffs.len(), n_bits);

    // extract cs from x_var
    let cs = x_var.cs();

    // --- 1) Extract witness value of x
    let x_val: F = x_var.value()?;
    let mut bits: Vec<bool> = x_val.into_bigint().to_bits_le();
    if bits.len() < n_bits {
        bits.resize(n_bits, false);
    } else {
        bits.truncate(n_bits);
    }

    // --- 2) Allocate bit witnesses
    let mut bit_vars = Vec::with_capacity(n_bits);
    for b in bits.iter() {
        let b_var = Boolean::<F>::new_witness(cs.clone(), || Ok(*b))?;
        bit_vars.push(b_var);
    }

    // --- 3) Reconstruct x from its bit decomposition
    let mut sum = FpVar::<F>::zero();
    for (coef, b) in pow2_coeffs.iter().zip(bit_vars.iter()) {
        // Boolean is stored as FpVar so it doesn't add any constraints
        let b_as_fp: FpVar<F> = b.clone().into();
        sum += coef.clone() * b_as_fp;
    }

    // --- 4) Enforce that x == sum
    x_var.enforce_equal(&sum)?;

    Ok(())
}

/// Enforce x < limit by decomposing product = (limit - x) * limit
/// into bits represented as FpVar and checking:
///    2^i * bit_i == product
///
/// - pow2_coeffs: expected length = number of bits (e.g., 2*q),
///   precomputed [1, 2, 4, ...] as FpVar constants.
pub fn enforce_x_less_than_limit<F: PrimeField>(
    x_var: &FpVar<F>,
    limit: &FpVar<F>,
    pow2_coeffs: &[FpVar<F>],
) -> Result<(), SynthesisError> {
    // --- 1) Witness value of product
    let x_val: F = x_var.value()?;
    let limit_val: F = limit.value()?;
    let product_val: F = (limit_val - x_val) * limit_val;

    // Convert product witness into bits
    let mut bits: Vec<bool> = product_val.into_bigint().to_bits_le();
    let n_bits = pow2_coeffs.len();
    if bits.len() < n_bits {
        bits.resize(n_bits, false);
    } else {
        bits.truncate(n_bits);
    }

    // --- 2) Allocate each bit as FpVar witness
    let cs = x_var.cs();
    let mut bit_vars: Vec<Boolean<F>> = Vec::with_capacity(n_bits);
    for &b in bits.iter() {
        let b_var = Boolean::<F>::new_witness(cs.clone(), || Ok(b))?;
        // (&b_var * &b_var - &b_var).enforce_equal(&FpVar::zero())?;
        bit_vars.push(b_var);
    }

    // --- 3) Build the sum inside the circuit
    let mut sum = FpVar::<F>::zero();
    for (coef, b) in pow2_coeffs.iter().zip(bit_vars.iter()) {
        // coef * b as a linear combination
        // Convert boolean to FpVar without new witness
        let b_as_fp: FpVar<F> = b.clone().into();
        sum += coef.clone() * b_as_fp;
    }

    // --- 4) Compute product in-circuit and enforce consistency
    let prod_var = (limit - x_var) * limit;
    prod_var.enforce_equal(&sum)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr;
    use ark_ff::AdditiveGroup;
    use ark_r1cs_std::fields::fp::FpVar;
    use ark_relations::r1cs::{ConstraintLayer, ConstraintSystem};
    use tracing_subscriber::layer::SubscriberExt;
    use tracing_subscriber::util::SubscriberInitExt;
    use ark_ff::{PrimeField};

    #[test]
    fn test_field_to_u128_conversion() {
        // Example values to test
        let test_values: [u128; 5] = [0, 1, 2, 15, 12345678901234567890];

        for &v in test_values.iter() {
            // Convert u128 -> Fr
            let f_val = Fr::from(v);

            // Convert Fr -> u128 using your code
            let big_int = f_val.into_bigint();
            let converted: u128 = big_int
                .to_bytes_le()
                .iter()
                .rev()
                .fold(0u128, |acc, &b| (acc << 8) | b as u128);

            // Assert equality
            assert_eq!(v, converted, "Failed for value {}", v);
        }
    }


    #[test]
    fn test_smallness_check_under_limit() {
        // Enable constraint logging (optional, for debugging)
        let _ = tracing_subscriber::registry()
            .with(ConstraintLayer::default())
            .try_init();

        // Build a new constraint system
        let cs = ConstraintSystem::<Fr>::new_ref();

        // Define limit = 20
        let limit: usize = 128;

        // Choose x < limit
        let x_val: u64 = (1 << 19) + 123; // ~ half of limit
        let x_f = Fr::from(x_val);

        // Allocate limit and x
        let x_var = FpVar::<Fr>::new_witness(cs.clone(), || Ok(x_f)).unwrap();

        // Precompute powers of two up to 2*20 bits = 40
        let n_bits = limit;
        let mut pow2_coeffs = Vec::with_capacity(n_bits);
        let mut cur = Fr::from(1u64);
        for _ in 0..n_bits {
            pow2_coeffs.push(FpVar::<Fr>::constant(cur.clone()));
            cur.double_in_place();
        }

        // Enforce the smallness condition
        enforce_x_less_than_limit_power_of_two(&x_var, limit, &pow2_coeffs).unwrap();

        // Number of constraints
        let num_constraints = cs.num_constraints();
        // Number of instance (public input) variables
        let num_instance = cs.num_instance_variables();
        // Get the witness assignment length directly
        let witness_size = cs.borrow().unwrap().witness_assignment.len();

        println!("constraints: {}", num_constraints);
        println!("public inputs: {}", num_instance);
        println!("witness size: {}", witness_size);

        // Check satisfiability
        assert!(cs.is_satisfied().unwrap());
    }

    #[test]
    fn test_enforce_x_less_than_limit() {
        // Setup tracing for debugging (optional)
        let _ = tracing_subscriber::registry().with(ConstraintLayer::default()).try_init();

        // Create a constraint system
        let cs = ConstraintSystem::<Fr>::new_ref();

        // Define a small limit (say limit = 8 = 2^3)
        let limit_val = Fr::from(8u64);
        let limit_var = FpVar::<Fr>::new_constant(cs.clone(), limit_val).unwrap();

        // Precompute pow2 coeffs [1,2,4,8,...]
        let pow2_coeffs: Vec<FpVar<Fr>> = (0..6)
            .map(|i| {
                let val = Fr::from(1u64 << i);
                FpVar::<Fr>::new_constant(cs.clone(), val).unwrap()
            })
            .collect();

        // Case 1: valid x < limit
        let x_val = Fr::from(5u64);
        let x_var = FpVar::<Fr>::new_witness(cs.clone(), || Ok(x_val)).unwrap();
        enforce_x_less_than_limit(&x_var, &limit_var, &pow2_coeffs).unwrap();

        assert!(cs.is_satisfied().unwrap());

        // Case 2: invalid x >= limit (x=9, limit=8)
        let cs2 = ConstraintSystem::<Fr>::new_ref();
        let limit_var2 = FpVar::<Fr>::new_constant(cs2.clone(), limit_val).unwrap();
        let pow2_coeffs2: Vec<FpVar<Fr>> = (0..4)
            .map(|i| {
                let val = Fr::from(1u64 << i);
                FpVar::<Fr>::new_constant(cs2.clone(), val).unwrap()
            })
            .collect();
        let x_val2 = Fr::from(9u64);
        let x_var2 = FpVar::<Fr>::new_witness(cs2.clone(), || Ok(x_val2)).unwrap();

        enforce_x_less_than_limit(&x_var2, &limit_var2, &pow2_coeffs2).unwrap();

        assert!(!cs2.is_satisfied().unwrap());
    }
}
