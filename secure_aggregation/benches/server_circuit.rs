use kzh_fold::constant_for_curves::{ScalarField as F, C2, E, G1, G2};
use kzh_fold::kzh::kzh3::{KZH3, KZH3SRS};
use kzh_fold::kzh::KZH;
use kzh_fold::kzh_fold::kzh3_fold::Accumulator3 as Accumulator;
use kzh_fold::nexus_spartan::commitment_traits::ToAffine;
use kzh_fold::nexus_spartan::committed_relaxed_snark::CRSNARKKey;
use kzh_fold::nexus_spartan::crr1cs::{is_sat, produce_synthetic_crr1cs, CRR1CSInstance, CRR1CSShape, CRR1CSWitness};
use kzh_fold::nexus_spartan::crr1csproof::CRR1CSProof;
use kzh_fold::nexus_spartan::matrix_evaluation_accumulation::verifier_circuit::{MatrixEvaluationAccVerifier, MatrixEvaluationAccVerifierVar};
use kzh_fold::nexus_spartan::partial_verifier::partial_verifier::SpartanPartialVerifier;
use kzh_fold::nexus_spartan::partial_verifier::partial_verifier_var::SpartanPartialVerifierVar;
use kzh_fold::nova::cycle_fold::coprocessor::setup_shape;
use kzh_fold::transcript::transcript::Transcript;
use kzh_fold::transcript::transcript_var::TranscriptVar;
use ark_ec::pairing::Pairing;
use ark_r1cs_std::alloc::{AllocVar, AllocationMode};
use ark_relations::r1cs::{ConstraintSystem, SynthesisMode};
use ark_serialize::CanonicalSerialize;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use kzh_fold::kzh3_verifier_circuit::prover::KZH3VerifierCircuitProver;
use kzh_fold::kzh3_verifier_circuit::verifier_circuit::KZH3VerifierVar;
use kzh_fold::nexus_spartan::conversion::convert_crr1cs;
use secure_aggregation::server_circuit::KZH3ServerCircuitVar;

fn bench(c: &mut Criterion) {
    // Directly map poseidon_num to (num_vars, num_inputs)
    let (num_vars, num_inputs) = (2 * 131072, 11);

    let (pcs_srs, spartan_shape, spartan_instance, spartan_proof, rx, ry) = {
        let num_cons = num_vars;

        // this generates a new instance/witness for spartan as well as PCS parameters
        let (spartan_shape, spartan_instance, spartan_witness, spartan_key) = produce_synthetic_crr1cs::<E, KZH3<E>>(num_cons, num_vars, num_inputs);

        assert!(is_sat(&spartan_shape, &spartan_instance, &spartan_witness, &spartan_key.gens_r1cs_sat).unwrap());

        let pcs_srs = spartan_key.gens_r1cs_sat.clone();

        let mut prover_transcript = Transcript::new(b"example");

        // Get `proof_i` and random evaluation point (r_x, r_y)
        let (spartan_proof, rx, ry) = CRR1CSProof::prove(
            &spartan_shape,
            &spartan_instance,
            spartan_witness,
            &spartan_key.gens_r1cs_sat,
            &mut prover_transcript,
        );

        (pcs_srs, spartan_shape, spartan_instance, spartan_proof, rx, ry)
    };

    // fresh transcripts to be used by the prover and verifier
    let mut prover_transcript = Transcript::new(b"example");
    let verifier_transcript_clone = prover_transcript.clone();
    let cs = ConstraintSystem::<F>::new_ref();

    let partial_verifier_var = {
        let mut verifier_transcript = prover_transcript.clone();
        // Get A(r_x, r_y), B(r_x, r_y), C(r_x, r_y)
        let current_A_B_C_evaluations = spartan_shape.inst.inst.evaluate(&rx, &ry);

        let partial_verifier = SpartanPartialVerifier::initialise(
            &spartan_proof,
            spartan_shape.get_num_vars(),
            spartan_shape.get_num_cons(),
            (spartan_instance.input.assignment, {
                let com_w: <E as Pairing>::G1Affine = spartan_instance.comm_W.clone().to_affine();
                com_w
            }),
            &current_A_B_C_evaluations,
            &mut prover_transcript,
        );

        partial_verifier.verify(&mut verifier_transcript);

        let partial_verifier_var = SpartanPartialVerifierVar::new_variable(
            cs.clone(),
            || Ok(partial_verifier.clone()),
            AllocationMode::Input,
        ).unwrap();

        partial_verifier_var
    };

    let acc_verifier_var = {
        let acc_srs = Accumulator::setup(pcs_srs.clone(), &mut thread_rng());

        // Get the KZH opening proof from the Spartan proof
        let opening_proof = spartan_proof.proof_eval_vars_at_ry.clone();

        // Commitment to witness polynomial
        let commitment_w = spartan_instance.comm_W.clone();

        let input = &ry[1..];

        // Sanity check: verify the opening proof
        KZH3::verify(
            &pcs_srs,
            &input,
            &spartan_proof.eval_vars_at_ry,
            &commitment_w,
            &opening_proof,
        );

        // Get accumulator from the opening proof
        let acc_instance = Accumulator::proof_to_accumulator_instance(
            &acc_srs,
            &input,
            &spartan_proof.eval_vars_at_ry,
            &commitment_w,
            &opening_proof,
        );

        let acc_witness = Accumulator::proof_to_accumulator_witness(
            &acc_srs,
            commitment_w,
            opening_proof,
            &input,
        );

        let current_acc = Accumulator::new(&acc_instance, &acc_witness);

        // println!("proof size: {}", proof.compressed_size());
        println!("acc size: {}", current_acc.compressed_size());

        // Check that the accumulator is valid
        Accumulator::decide(
            &acc_srs,
            &current_acc,
        );

        // use a random accumulator as the running one
        let running_acc = Accumulator::rand(&acc_srs);

        // the shape of the R1CS instance
        let ova_shape = setup_shape::<G1, G2>().unwrap();

        // get trivial running instance
        let (ova_running_instance, ova_running_witness) = KZH3VerifierCircuitProver::<G1, G2, C2, E, F>::get_trivial_cycle_fold_running_instance_witness(&ova_shape);

        // get commitment_pp
        let ova_commitment_pp = KZH3VerifierCircuitProver::<G1, G2, C2, E, F>::get_commitment_pp(&ova_shape);

        let kzh_acc_verifier_prover: KZH3VerifierCircuitProver<G1, G2, C2, E, F> = KZH3VerifierCircuitProver::new(
            &acc_srs,
            ova_commitment_pp,
            running_acc,
            current_acc.clone(),
            ova_running_instance,
            ova_running_witness,
            prover_transcript,
        );

        // assert it's formated correctly
        kzh_acc_verifier_prover.is_satisfied();

        let acc_verifier_var = KZH3VerifierVar::<G1, G2, C2>::new::<E>(cs.clone(), kzh_acc_verifier_prover);

        acc_verifier_var
    };

    let matrix_evaluation_verifier_var = {
        let matrix_eval_acc_verifier = MatrixEvaluationAccVerifier::random_from_eval_point(
            &spartan_shape,
            rx.clone(),
            ry.clone(),
            &mut thread_rng(),
        );

        let matrix_evaluation_verifier_var = MatrixEvaluationAccVerifierVar::new_variable(
            cs.clone(),
            || Ok(matrix_eval_acc_verifier.clone()),
            AllocationMode::Input,
        ).unwrap();

        matrix_evaluation_verifier_var
    };

    // construct the augmented circuit
    let augmented_circuit = KZH3ServerCircuitVar {
        spartan_partial_verifier: partial_verifier_var,
        kzh_acc_verifier: acc_verifier_var,
        matrix_evaluation_verifier: matrix_evaluation_verifier_var,
    };

    let mut transcript_var = TranscriptVar::from_transcript(cs.clone(), verifier_transcript_clone);

    // run the verification function on augmented circuit
    let _ = augmented_circuit.verify::<E>(&pcs_srs, cs.clone(), &mut transcript_var);

    println!("augmented circuit constraints: {}", cs.num_constraints());

    // Set the mode to Prove before we convert it for spartan
    cs.set_mode(SynthesisMode::Prove { construct_matrices: true });
    cs.finalize();

    // convert to the corresponding Spartan types
    let shape = CRR1CSShape::<F>::convert::<G1>(cs.clone());

    // get the number the minimum size we need for committing to the constraint system
    let min_num_vars = CRSNARKKey::<E, KZH3<E>>::get_min_num_vars(shape.get_num_cons(), shape.get_num_vars(), shape.get_num_inputs());
    let srs: KZH3SRS<E> = KZH3::setup(min_num_vars + 1, &mut thread_rng());


    let bench_name = "committing to the polynomial".to_string();
    c.bench_function(&bench_name, |b| {
        b.iter(|| {
            let (_, _): (CRR1CSInstance<E, KZH3<E>>, CRR1CSWitness<E, KZH3<E>>) = convert_crr1cs(cs.clone(), &srs);
        })
    });

    let (instance, witness): (CRR1CSInstance<E, KZH3<E>>, CRR1CSWitness<E, KZH3<E>>) = convert_crr1cs(cs.clone(), &srs);

    let bench_name = "run spartan proof".to_string();
    c.bench_function(&bench_name, |b| {
        b.iter(|| {
            let mut new_prover_transcript = Transcript::new(b"example");
            let (_, rx, ry) = CRR1CSProof::prove(
                &shape,
                &instance,
                witness.clone(),
                &srs,
                &mut new_prover_transcript,
            );

            // evaluate matrices A B C
            let _ = shape.inst.inst.evaluate(&rx, &ry);
        })
    });

}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = server_circuit;
    config = custom_criterion_config();
    targets = bench
}

criterion_main!(server_circuit);