use ark_ff::{AdditiveGroup, One};
use ark_std::UniformRand;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use kzh_fold::constant_for_curves::{ScalarField, E};
use kzh_fold::kzh::KZH;
use kzh_fold::kzh::kzh4::KZH4;
use kzh_fold::nexus_spartan::sumcheck::SumcheckInstanceProof;
use kzh_fold::polynomial::eq_poly::eq_poly::EqPolynomial;
use kzh_fold::polynomial::multilinear_poly::multilinear_poly::MultilinearPolynomial;
use kzh_fold::transcript::transcript::Transcript;

fn bench(c: &mut Criterion) {
    let num_variables = vec![10, 17, 23];
    for n in num_variables {
        let srs = KZH4::<E>::setup(n, &mut thread_rng());

        // generate random binary polynomial
        let mut poly: MultilinearPolynomial<ScalarField> = MultilinearPolynomial::random_binary(n, &mut thread_rng());

        // commit to the polynomial
        let (com, aux) = KZH4::commit(&srs, &poly, &mut thread_rng());

        // get a random challenge r from the verifier, can be derived through Fiat-Shamir too
        let r: Vec<ScalarField> = (0..n)
            .map(|_| ScalarField::rand(&mut thread_rng()))
            .collect();

        let eq = EqPolynomial::new(r);
        let eq_poly = MultilinearPolynomial::new(eq.evals());

        let bench_name = format!("server time, n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let eq = EqPolynomial::new(r);
                let eq_poly = MultilinearPolynomial::new(eq.evals());

                // run the sumcheck, first server cost
                let (_, eval_point, _) =  SumcheckInstanceProof::prove_cubic_three_terms(
                    &ScalarField::ZERO,
                    n,
                    &mut poly.clone(),
                    &mut poly.clone(),
                    &mut eq_poly.clone(),
                    |a: &ScalarField, b: &ScalarField, c: &ScalarField| -> ScalarField { a.clone() * (ScalarField::one() - b.clone()) * c.clone() },
                    &mut Transcript::new(b"new transcript")
                );

                // second server cost
                let _ = KZH4::open(&srs, eval_point.as_slice(), &com, &aux, &poly, &mut thread_rng());
            })
        });

        // run the sumcheck, first server cost
        let (sumcheck_proof, eval_point, eval_res) =  SumcheckInstanceProof::prove_cubic_three_terms(
            &ScalarField::ZERO,
            n,
            &mut poly.clone(),
            &mut poly.clone(),
            &mut eq_poly.clone(),
            |a: &ScalarField, b: &ScalarField, c: &ScalarField| -> ScalarField { a.clone() * (ScalarField::one() - b.clone()) * c.clone() },
            &mut Transcript::new(b"new transcript")
        );

        // second server cost
        let open = KZH4::open(&srs, eval_point.as_slice(), &com, &aux, &poly, &mut thread_rng());

        let bench_name = format!("client time, n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                // verify sumcheck
                sumcheck_proof.verify::<E>(ScalarField::ZERO, n, 3, &mut Transcript::new(b"new transcript")).expect("TODO: panic message");

                // first client cost
                assert_eq!(eq.evaluate(eval_point.as_slice()), eval_res[2].clone());

                // second client cost
                KZH4::verify(&srs, eval_point.as_slice(), &eval_res[0].clone(), &com, &open)
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = uniqueness_proof;
    config = custom_criterion_config();
    targets = bench
}

criterion_main!(uniqueness_proof);