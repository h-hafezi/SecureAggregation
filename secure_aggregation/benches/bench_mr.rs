use criterion::{criterion_group, criterion_main, Criterion};
use rand::thread_rng;
use kzh_fold::constant_for_curves::{G1Affine, E};
use kzh_fold::kzh::KZH;
use kzh_fold::kzh::kzh4::{KZH4, KZH4SRS};
use secure_aggregation::merkle_tree::mr::MerkleTree;

fn bench(c: &mut Criterion) {
    let num_variables = vec![10, 20, 21, 30];
    for n in num_variables {
        // get srs
        let srs: KZH4SRS<E> = KZH4::setup(n, &mut thread_rng());


        let bench_name = format!("Bench KZH4 SRS n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                // bench srs
                let _: KZH4SRS<E> = KZH4::setup(n, &mut thread_rng());
            })
        });

        // generate the MR
        let flatten_array: Vec<G1Affine> = srs.H_xyzt.into_iter().collect();
        assert_eq!(flatten_array.len(), 1 << n);

        let bench_name = format!("Build tree parallel with 4 threads, n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = MerkleTree::<E>::build_parallel(flatten_array.clone(), 4);
            })
        });


        let bench_name = format!("Build tree parallel with 8 threads, n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = MerkleTree::<E>::build_parallel(flatten_array.clone(), 8);
            })
        });

        let bench_name = format!("Build tree parallel with 16 threads, n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = MerkleTree::<E>::build_parallel(flatten_array.clone(), 16);
            })
        });

        let bench_name = format!("Build tree parallel with 32 threads, n={}", n);
        c.bench_function(&bench_name, |b| {
            b.iter(|| {
                let _ = MerkleTree::<E>::build_parallel(flatten_array.clone(), 32);
            })
        });
    }
}

fn custom_criterion_config() -> Criterion {
    Criterion::default().sample_size(10)
}

// Benchmark group setup
criterion_group! {
    name = bench_mr;
    config = custom_criterion_config();
    targets = bench
}

criterion_main!(bench_mr);