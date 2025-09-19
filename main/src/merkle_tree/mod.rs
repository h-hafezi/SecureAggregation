use crate::gadgets::non_native::non_native_affine_var::NonNativeAffineVar;
use crate::hash::poseidon::PoseidonHashVar;
use ark_crypto_primitives::sponge::constraints::{AbsorbGadget, CryptographicSpongeVar};
use ark_crypto_primitives::sponge::Absorb;
use ark_ec::CurveConfig;
use ark_ff::PrimeField;
use ark_r1cs_std::alloc::AllocVar;
use ark_r1cs_std::fields::fp::FpVar;
use ark_r1cs_std::R1CSVar;
use ark_relations::r1cs::SynthesisError;
use ark_std::UniformRand;

pub fn merkle_tree<G1>(
    poseidon_var: &mut PoseidonHashVar<G1::ScalarField>,
    leaf: &NonNativeAffineVar<G1>,
    depth: usize,
) -> Result<FpVar<G1::ScalarField>, SynthesisError>
where
    G1: ark_ec::short_weierstrass::SWCurveConfig,
    G1::BaseField: PrimeField, <G1 as CurveConfig>::ScalarField: Absorb
{
    let cs = poseidon_var.sponge.cs();

    // 1. Convert leaf into FpVars
    let leaf_elems: Vec<FpVar<G1::ScalarField>> = leaf.to_sponge_field_elements()?;

    // 2. Absorb leaf into the Poseidon sponge
    poseidon_var.update_sponge(leaf_elems);

    // 3. Initial sponge output
    let mut output = poseidon_var.output();

    // 4. Loop for depth
    for _ in 0..depth {
        // Generate a witness FpVar with a "random" value (deterministic in circuit)
        let rand_var = FpVar::<G1::ScalarField>::new_witness(cs.clone(), || {
            Ok(G1::ScalarField::rand(&mut ark_std::test_rng()))
        })?;

        // Absorb into sponge
        poseidon_var.update_sponge(vec![rand_var]);

        // Optionally, you can refresh output in each round if desired:
        output = poseidon_var.output();
    }

    // 5. Return final output
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_r1cs_std::R1CSVar;
    use ark_relations::r1cs::ConstraintSystem;
    use ark_std::test_rng;

    use crate::constant_for_curves::{G1Projective, ScalarField, G1};

    #[test]
    fn merkle_tree_test() {
        // === Native side ===
        // Pick a random affine point as the "leaf"
        let mut rng = test_rng();
        let point_native: G1Projective = G1Projective::rand(&mut rng);
        let affine_native = point_native.into_affine();

        // === Circuit side ===
        let cs = ConstraintSystem::<ScalarField>::new_ref();

        // Allocate the leaf in the circuit
        let leaf_var = NonNativeAffineVar::<G1>::new_witness(
            cs.clone(),
            || Ok(point_native),
        ).unwrap();

        // Create PoseidonVar
        let mut poseidon_var = PoseidonHashVar::<ScalarField>::new(cs.clone());

        // Run the Merkle tree gadget
        let depth = 20;
        let out_var = merkle_tree::<G1>(&mut poseidon_var, &leaf_var, depth).unwrap();

        // Check constraint system health
        assert!(cs.is_satisfied().unwrap());
        println!("num_constraints: {}", cs.num_constraints());

        // Get the circuit-evaluated value
        let out_value = out_var.value().unwrap();
    }
}

