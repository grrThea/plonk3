use p3_baby_bear::{BabyBear, DiffusionMatrixBabybear};
use core::mem::{size_of};
use std::clone;

use p3_air::{Air, AirBuilder, BaseAir};
use p3_challenger::DuplexChallenger;
use p3_matrix::dense::RowMajorMatrix;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::Field;
use p3_matrix::{Matrix, MatrixRowSlices};
use p3_poseidon2::Poseidon2;
use p3_uni_stark::{prove, verify, StarkConfig, StarkGenericConfig, Val, VerificationError};
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};
use tracing_forest::util::LevelFilter;
use tracing_forest::ForestLayer;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_commit::ExtensionMmcs;
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_util::log2_ceil_usize;


pub struct FibonacciAir {
    a: u64,
    b: u64,
    n: usize,
    x: usize,
}

impl Default for FibonacciAir {
    fn default() -> Self {
        Self {
            a: 0,
            b: 0,
            n: 0,
            x: 0,
        }
    }
}

#[derive(Clone)]
struct FibonacciCol {
    pre: u64,
    cur: u64,
    sum: u64,
    is_last: u64,
}

impl Default for FibonacciCol {
    fn default() -> Self {
        Self {
            pre: 0,
            cur: 0,
            sum: 0,
            is_last: 0,
        }
    }
}

pub(crate) const NUM_FIB_COLS: usize = size_of::<FibonacciCol>();


impl FibonacciAir {
    pub fn random_valid_trace<F: Field>(&self) -> RowMajorMatrix<F>
    where
        Standard: Distribution<F>,
    {
        let mut trace_values = vec![F::default(); 4 * 2];
        trace_values[0] = F::from_wrapped_u64(self.a);
        trace_values[1] = F::from_wrapped_u64(self.b);
        for i in 2..self.n as usize {
            trace_values[i] = trace_values[i-1] + trace_values[i-2];
        }
        println!("trace_values =>{:?}", trace_values);
       

        RowMajorMatrix::new(trace_values, 2)
    }
}

impl<F> BaseAir<F> for FibonacciAir {
    fn width(&self) -> usize {
        2
    }
}

impl<AB: AirBuilder> Air<AB> for FibonacciAir {
    fn eval(&self, builder: &mut AB) {
        let main: <AB as AirBuilder>::M = builder.main();
        let local = main.row_slice(0);
        let next = main.row_slice(1);
        builder.when_first_row().assert_one(local[0]);
        builder.when_first_row().assert_one(local[1]);
        builder.when_transition().assert_eq(local[0]+local[1], next[0]);
        builder.when_transition().assert_eq(local[1]+next[0], next[1]);

        // builder.when_last_row().assert_eq(local[1]+local[0], self.x);
    }
}


fn main() -> Result<(), VerificationError> {
  let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();

    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .init();

        type Val = BabyBear;
    type Challenge = BinomialExtensionField<Val, 4>;

    type Perm = Poseidon2<Val, DiffusionMatrixBabybear, 16, 7>;
    let perm = Perm::new_from_rng(8, 22, DiffusionMatrixBabybear, &mut thread_rng());

    type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
    let hash = MyHash::new(perm.clone());

    type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
    let compress = MyCompress::new(perm.clone());

    type ValMmcs = FieldMerkleTreeMmcs<
        <Val as Field>::Packing,
        <Val as Field>::Packing,
        MyHash,
        MyCompress,
        8,
    >;
    let val_mmcs = ValMmcs::new(hash, compress);

    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    type Dft = Radix2DitParallel;
    let dft = Dft {};

    type Challenger = DuplexChallenger<Val, Perm, 16>;

    let a = 1;
    let b = 1;
    let n = 8;
    let x = 21; 
    let air = FibonacciAir { a, b, n, x };

    let trace = air.random_valid_trace();
    println!("trace: {:?}", trace);

    let fri_config = FriConfig {
        log_blowup: 1,
        num_queries: 100,
        proof_of_work_bits: 16,
        mmcs: challenge_mmcs,
    };
    type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
    let pcs = Pcs::new(log2_ceil_usize(trace.height()), dft, val_mmcs, fri_config);

    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    let config = MyConfig::new(pcs);

    let mut challenger = Challenger::new(perm.clone());

    let proof = prove::<MyConfig, _>(&config, &FibonacciAir::default(), &mut challenger, trace);

    let mut challenger = Challenger::new(perm);
    verify(&config, &FibonacciAir::default(), &mut challenger, &proof)
}

