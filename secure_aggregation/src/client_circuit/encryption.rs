use ark_ff::{BigInteger, PrimeField, Zero};
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::eq::EqGadget;
use ark_r1cs_std::fields::FieldVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::{ConstraintSystemRef, SynthesisError};
use num::One;
use num_bigint::BigUint;
use rand::Rng;
use crate::client_circuit::smallness::{enforce_x_less_than_limit, enforce_x_less_than_limit_power_of_two};

/// The secure_aggregation composite circuit function.
///
/// - `r1` and `r2` are FpVar witnesses, and `q`, `ell` are usize (r1, r2 random moduli q).
/// - `v` is a witness vector of length N (values mod ell).
/// - `u` is a public vector of length N (constants allocated outside, passed as FpVar::constant).
/// - `res1` and `res2` are witness vectors of length N (values in F_q).
/// - `pow2_q` and `pow2_ell` are coefficient arrays prepared for the smallness checks.
///
pub fn AHE_relation<F: PrimeField>(
    r1: &FpVar<F>,
    r2: &FpVar<F>,
    q: usize, // log of the modulo is in fact provided
    ell: &FpVar<F>,
    v: &[FpVar<F>], // witness which is mod ell
    u: &[FpVar<F>], // public constants
    res: &[FpVar<F>], // witness, elements in F_q
    pow2_q: &[FpVar<F>], // coefficients for smallness relative to q of length 2 * q
    pow2_ell: &[FpVar<F>], // coefficients for smallness relative to ell of length 2 * log(ell)
) -> Result<(), SynthesisError> {
    let cs = r1.cs();
    let q_var = FpVar::Constant(F::from_be_bytes_mod_order((BigUint::one() << q).to_bytes_be().as_slice()));
    let n = v.len();

    // --- 1) Ensure all elements in v are less than ell
    for val_v in v {
        enforce_x_less_than_limit(val_v, ell, pow2_ell)?;
    }

    // --- 2) Ensure res1 and res2 are all less than q
    for val_res1 in res {
        enforce_x_less_than_limit_power_of_two(val_res1, q, pow2_q)?;
    }

    // --- 3) Compute powers of r1 modulo q and enforce correctness
    let mut r1_powers: Vec<FpVar<F>> = Vec::with_capacity(n);
    // Base case: r1^0 = 1
    let mut cur_r1_power_var = FpVar::<F>::constant(F::one());
    r1_powers.push(cur_r1_power_var.clone());


    for i in 1..n {
        let prev_r1_power_var = r1_powers[i - 1].clone();

        // Witness computation (outside the circuit)
        let prod_val = prev_r1_power_var.value()? * r1.value()?;
        let prod_bigint = prod_val.into_bigint();
        let prod = BigUint::from_bytes_le(&prod_bigint.to_bytes_le());
        let q_bigint = BigUint::one() << q;


        // compute quotient and remainder
        let quotient = &prod / &q_bigint;
        let remainder = &prod % &q_bigint;

        // convert back into ArkWorks BigInt
        let m = F::from_le_bytes_mod_order(quotient.to_bytes_le().as_ref());
        let r1_power_i = F::from_le_bytes_mod_order(remainder.to_bytes_le().as_ref());

        // Allocate quotient (m_i) and new power (r1^i) as witnesses
        let m_var = FpVar::<F>::new_witness(cs.clone(), || Ok(m))?;
        let r1_power_i_var = FpVar::<F>::new_witness(cs.clone(), || Ok(r1_power_i))?;

        // Enforce the remainder is less than q
        enforce_x_less_than_limit_power_of_two(&r1_power_i_var, q, pow2_q)?;

        // Enforce the modular relation: prev_power * r1 = q * m + r1^i
        let left = prev_r1_power_var * r1;
        let right = &q_var * m_var + r1_power_i_var.clone();
        left.enforce_equal(&right)?;

        // Add the new power to the vector
        r1_powers.push(r1_power_i_var);
    }

    // --- 3) Compute powers of r2 modulo q and enforce correctness
    let mut r2_powers: Vec<FpVar<F>> = Vec::with_capacity(n);
    // Base case: r2^0 = 1
    let mut cur_r2_power_var = FpVar::<F>::constant(F::one());
    r2_powers.push(cur_r2_power_var.clone());

    for i in 1..n {
        let prev_r2_power_var = r2_powers[i - 1].clone();

        // Witness computation (outside the circuit)
        let prod_val = prev_r2_power_var.value()? * r2.value()?;
        let prod_bigint = prod_val.into_bigint();
        let prod = BigUint::from_bytes_le(&prod_bigint.to_bytes_le());
        let q_bigint = BigUint::one() << q;

        // compute quotient and remainder
        let quotient = &prod / &q_bigint;
        let remainder = &prod % &q_bigint;

        // convert back into ArkWorks BigInt
        let m = F::from_le_bytes_mod_order(quotient.to_bytes_le().as_ref());
        let r2_power_i = F::from_le_bytes_mod_order(remainder.to_bytes_le().as_ref());

        // Allocate quotient (m_i) and new power (r2^i) as witnesses
        let m_var = FpVar::<F>::new_witness(cs.clone(), || Ok(m))?;
        let r2_power_i_var = FpVar::<F>::new_witness(cs.clone(), || Ok(r2_power_i))?;

        // Enforce the remainder is less than q
        enforce_x_less_than_limit_power_of_two(&r2_power_i_var, q, pow2_q)?;

        // Enforce the modular relation: prev_power * r2 = q * m + r2^i
        let left = prev_r2_power_var * r2;
        let right = &q_var * m_var + r2_power_i_var.clone();
        left.enforce_equal(&right)?;

        // Add the new power to the vector
        r2_powers.push(r2_power_i_var);
    }

    // --- 4.2) Compute running sums V'_n, U'_n, R2_n for r2
    let mut v_r1 = FpVar::<F>::constant(F::zero());
    let mut u_r1 = FpVar::<F>::constant(F::zero());
    let mut res_r1 = FpVar::<F>::constant(F::zero());

    for i in 0..n {
        let r_i_power = r1_powers[i].clone();

        // V'_i = V'_{i-1} + v[i] * r2^i
        v_r1 += v[i].clone() * r_i_power.clone();

        // U'_i = U'_{i-1} + u[i] * r2^i
        u_r1 += u[i].clone() * r_i_power.clone();

        // R2_i = R2_{i-1} + res2[i] * r2^i
        res_r1 += res[i].clone() * r_i_power;
    }

    // --- 4.2) Compute running sums V'_n, U'_n, R2_n for r2
    let mut v_r2 = FpVar::<F>::constant(F::zero());
    let mut u_r2 = FpVar::<F>::constant(F::zero());
    let mut res_r2 = FpVar::<F>::constant(F::zero());

    for i in 0..n {
        let r_i_power = r2_powers[i].clone();

        // V'_i = V'_{i-1} + v[i] * r2^i
        v_r2 += v[i].clone() * r_i_power.clone();

        // U'_i = U'_{i-1} + u[i] * r2^i
        u_r2 += u[i].clone() * r_i_power.clone();

        // R2_i = R2_{i-1} + res2[i] * r2^i
        res_r2 += res[i].clone() * r_i_power;
    }

    Ok(())
}

pub fn get_AHE_params<F:PrimeField>(cs: ConstraintSystemRef<F>) -> (
    FpVar<F>,
    FpVar<F>,
    usize,
    FpVar<F>,
    Vec<FpVar<F>>,
    Vec<FpVar<F>>,
    Vec<FpVar<F>>,
    Vec<FpVar<F>>,
    Vec<FpVar<F>>
) {
    const N: usize = 4096; // Vector length (Polynomial degree N-1)

    // Moduli values: Q = 2^96, ell = 2^20
    let ell_big: BigUint = BigUint::from(1u64) << 20;

    // Convert BigUint to Field elements
    let q = 96usize;
    let ell = F::from_le_bytes_mod_order(ell_big.to_bytes_le().as_ref());

    // Random challenge values (kept small for computation feasibility)
    let mut rng = rand::thread_rng();
    let r1_val: u64 = rng.gen_range(1..100);
    let r2_val: u64 = rng.gen_range(1..100);

    let r1 = F::from(r1_val);
    let r2 = F::from(r2_val);

    // Input vectors, must satisfy smallness v_i < ell
    let u_vals: Vec<u64> = (0..N).map(|_| rng.gen_range(1..1000)).collect();
    let v_vals: Vec<u64> = (0..N).map(|_| rng.gen_range(1..1000)).collect();

    // Compute the resulting polynomial
    let res_vals = vec![F::zero(); N];

    // --- 3. Allocate FpVar Witnesses/Constants ---
    let ell_var = FpVar::<F>::new_witness(cs.clone(), || Ok(ell)).unwrap();

    let r1_var = FpVar::<F>::new_witness(cs.clone(), || Ok(r1)).unwrap();
    let r2_var = FpVar::<F>::new_witness(cs.clone(), || Ok(r2)).unwrap();

    // u is public constant
    let u_vars: Vec<FpVar<F>> = u_vals.iter()
        .map(|&val| FpVar::<F>::constant(F::from(val)))
        .collect();

    // v, res are witnesses
    let v_vars: Vec<FpVar<F>> = v_vals.iter()
        .map(|&val| FpVar::<F>::new_witness(cs.clone(), || Ok(F::from(val))).unwrap())
        .collect();

    let res_vars: Vec<FpVar<F>> = res_vals.iter()
        .map(|&val| FpVar::<F>::new_witness(cs.clone(), || Ok(val)).unwrap())
        .collect();

    let pow2_q_vars = get_pow2_coeffs(cs.clone(), q);
    let pow2_ell_vars = get_pow2_coeffs(cs.clone(), ell_big.bits() as usize * 2);

    (
        r1_var, r2_var, q, ell_var,
        v_vars, u_vars, res_vars,
        pow2_q_vars, pow2_ell_vars
    )
}

/// Helper function to create powers of two constants for smallness checks.
pub fn get_pow2_coeffs<F: PrimeField>(cs: ConstraintSystemRef<F>, bit_len: usize) -> Vec<FpVar<F>> {
    let mut pow2_coeffs = Vec::with_capacity(bit_len);
    let mut cur = F::from(1u64);
    for _ in 0..bit_len {
        pow2_coeffs.push(FpVar::<F>::constant(cur));
        cur.double_in_place();
    }
    pow2_coeffs
}

/// Helper to compute polynomial evaluation $\sum a_i r^i$ in $F_P$
fn poly_eval<F: PrimeField>(a: &[u64], r: F) -> F {
    let mut sum = F::zero();
    let mut r_power = F::one();
    for &a_i in a.iter() {
        sum += F::from(a_i) * r_power;
        r_power *= r;
    }
    sum
}

/// Helper: smallest power of two >= n
fn next_power_of_two_usize(n: usize) -> usize {
    if n.is_power_of_two() {
        n
    } else {
        n.next_power_of_two()
    }
}

fn average_bit_length<F: PrimeField>(witness: &[F]) -> f64 {
    let mut total_bits = 0usize;
    for elem in witness {
        // Convert to BigUint
        let bigint = elem.into_bigint();
        let n: BigUint = bigint.into();

        let bits = if n.is_zero() {
            1
        } else {
            n.bits() as usize
        };
        total_bits += bits;
    }

    if witness.is_empty() {
        0.0
    } else {
        total_bits as f64 / witness.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use std::time::{Duration, Instant};
    use super::*;
    use ark_ec::CurveGroup;
    use ark_ec::short_weierstrass::{Affine, Projective};
    use ark_ff::Zero;
    use ark_relations::r1cs::{ConstraintSystem, ConstraintSystemRef};
    use ark_std::rand::Rng;
    use crate::constant_for_curves::{ScalarField, G1};
    use crate::pederson::PedersenCommitment;
    use crate::commitment::CommitmentScheme;

    #[test]
    fn test_large_scale_constraint_count() {
        // Build a new constraint system
        let cs = ConstraintSystem::<ScalarField>::new_ref();
        let mut rng = rand::thread_rng();

        // --- 1. Define Parameters and Witnesses ---
        const N: usize = 4096; // Vector length (Polynomial degree N-1)

        // Moduli values: Q = 2^96, ell = 2^20
        let q_big: BigUint = BigUint::from(1u64) << 96;
        let ell_big: BigUint = BigUint::from(1u64) << 20;
        println!("hhhhhhhhhhhhhhhhh {}", ell_big);

        // Convert BigUint to Field elements
        let q = 96usize;
        let ell = ScalarField::from_le_bytes_mod_order(ell_big.to_bytes_le().as_ref());

        // Random challenge values (kept small for computation feasibility)
        let r1_val: u64 = rng.gen_range(1..100);
        let r2_val: u64 = rng.gen_range(1..100);

        let r1 = ScalarField::from(r1_val);
        let r2 = ScalarField::from(r2_val);

        // Input vectors, must satisfy smallness v_i < ell
        let u_vals: Vec<u64> = (0..N).map(|_| rng.gen_range(1..1000)).collect();
        let v_vals: Vec<u64> = (0..N).map(|_| rng.gen_range(1..1000)).collect();

        // Compute the resulting polynomial
        let res_vals = vec![ScalarField::zero(); N];

        // --- 3. Allocate FpVar Witnesses/Constants ---
        let ell_var = FpVar::<ScalarField>::new_witness(cs.clone(), || Ok(ell)).unwrap();

        let r1_var = FpVar::<ScalarField>::new_witness(cs.clone(), || Ok(r1)).unwrap();
        let r2_var = FpVar::<ScalarField>::new_witness(cs.clone(), || Ok(r2)).unwrap();

        // u is public constant
        let u_vars: Vec<FpVar<ScalarField>> = u_vals.iter()
            .map(|&val| FpVar::<ScalarField>::constant(ScalarField::from(val)))
            .collect();

        // v, res are witnesses
        let v_vars: Vec<FpVar<ScalarField>> = v_vals.iter()
            .map(|&val| FpVar::<ScalarField>::new_witness(cs.clone(), || Ok(ScalarField::from(val))).unwrap())
            .collect();

        let res_vars: Vec<FpVar<ScalarField>> = res_vals.iter()
            .map(|&val| FpVar::<ScalarField>::new_witness(cs.clone(), || Ok(val)).unwrap())
            .collect();

        // Precomputed power of 2 coefficients
        let ell_bit_len = ell_big.bits() as usize;

        let pow2_q_vars = get_pow2_coeffs(cs.clone(), q);
        let pow2_ell_vars = get_pow2_coeffs(cs.clone(), ell_bit_len * 2);

        // --- 4. Run the Circuit ---
        let result = AHE_relation(
            &r1_var, &r2_var, q, &ell_var,
            &v_vars, &u_vars, &res_vars,
            &pow2_q_vars, &pow2_ell_vars
        );

        // --- 5. Output Constraint Count ---
        assert!(result.is_ok(), "Circuit initialization should succeed.");
        println!("Configuration: N={}, q=2^96, ell=2^20", N);
        println!("Final Constraint Count: {}", cs.num_constraints());

        // Final assertion to ensure the witness is valid (optional but good practice)
        assert!(cs.is_satisfied().unwrap(), "Constraint system must be satisfied by the generated witness.");

        // OR clone if you need ownership
        let witness = {
            let cs_borrow = cs.borrow();
            let cs_ref = cs_borrow.as_ref().unwrap();
            cs_ref.witness_assignment.clone()
        };

        let len = next_power_of_two_usize(witness.len());

        println!("witness length: {}", witness.len());
        println!("average bit length: {}", average_bit_length(&witness));

        // Setup witness and Pedersen commitment params
        let pp: Vec<Affine<G1>> = PedersenCommitment::<Projective<G1>>::setup(witness.len(), b"test", &());

        // 3) Benchmark
        let mut total_time = Duration::ZERO;
        for _ in 0..10 {
            let start = Instant::now();
            let _commitment: Projective<G1> = PedersenCommitment::commit(&pp, witness.as_slice());
            total_time += start.elapsed();
        }

        let avg_time = total_time / 10;
        println!(
            "Pedersen commit: witness.len = {}, repeats = {:?}, avg = {:?}",
            witness.len(),
            10,
            avg_time
        );
    }
}






