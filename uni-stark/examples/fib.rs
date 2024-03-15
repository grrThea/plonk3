use std::borrow::Borrow;

use p3_baby_bear::{BabyBear, DiffusionMatrixBabybear};

use p3_air::{Air, AirBuilder, BaseAir};
use p3_challenger::DuplexChallenger;
use p3_commit::testing::TrivialPcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::Field;
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::MatrixRowSlices;
use p3_poseidon2::Poseidon2;
use p3_uni_stark::{prove, verify, StarkConfig, StarkGenericConfig, Val, VerificationError};
use rand::distributions::{Distribution, Standard};
use rand::thread_rng;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};
use tracing_forest::util::LevelFilter;
use tracing_forest::ForestLayer;

pub struct FibonacciAir {
    a: u64,
    b: u64,
    n: usize,
    x: usize,
}

impl FibonacciAir {
    pub fn random_valid_trace<F: Field>(&self) -> RowMajorMatrix<F>
    where
        Standard: Distribution<F>,
    {
        let num_rows = self.n + 1;

        let mut trace_values = vec![F::default(); num_rows];
        trace_values[0] = F::from_wrapped_u64(self.a);
        trace_values[1] = F::from_wrapped_u64(self.b);
        for i in 2..self.n as usize {
            trace_values[i] = trace_values[i-1] + trace_values[i-2] 
        }
        trace_values[self.n] = F::from_canonical_usize(self.x);
        RowMajorMatrix::new(trace_values, 3)
    }
}

impl<F> BaseAir<F> for FibonacciAir {
    fn width(&self) -> usize {
        3
    }
}

impl<AB: AirBuilder> Air<AB> for FibonacciAir {
    fn eval(&self, builder: &mut AB) {
        let main: <AB as AirBuilder>::M = builder.main();
        let local = main.row_slice(0);

        let len= local.len();

        println!("len =>{:?}", len);


        for i in 2..self.n {
          let cur = local[i];
          let pre = local[i-1];
          let pre2 = local[i-2];
          builder.assert_eq(cur, pre + pre2)
        }

        let n: <AB as AirBuilder>::Var = local[len - 2];
        let x: <AB as AirBuilder>::Var = local[len - 1];
        builder.assert_eq(n,x);
    }
}

fn do_test<SC: StarkGenericConfig>(
    config: SC,
    a: u64,
    b: u64,
    n: usize,
    x: usize,
    challenger: SC::Challenger,
) -> Result<(), VerificationError>
where
    SC::Challenger: Clone,
    Standard: Distribution<Val<SC>>,
{
    let air = FibonacciAir { a, b, n, x };
    let trace = air.random_valid_trace();

    println!("trace => {:?}", trace);
    let mut p_challenger = challenger.clone();
    let proof = prove(&config, &air, &mut p_challenger, trace);

    let serialized_proof = postcard::to_allocvec(&proof).expect("unable to serialize proof");
    tracing::debug!("serialized_proof len: {} bytes", serialized_proof.len());

    let proof =
        postcard::from_bytes(&serialized_proof).expect("unable to deserialize proof");

    let mut v_challenger = challenger.clone();
    verify(&config, &air, &mut v_challenger, &proof)
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

    type Dft = Radix2DitParallel;
    let dft = Dft {};

    type Challenger = DuplexChallenger<Val, Perm, 16>;

    let log_n = 0;
    type Pcs = TrivialPcs<Val, Radix2DitParallel>;
    let pcs = Pcs::new(log_n, dft);

    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    let config = MyConfig::new(pcs);
    let a = 1;
    let b = 1;
    let n = 11;
    let x = 89; 
    // 1,1,2,3,5,8,13,21,34,55,89
    do_test(config, a, b, n, x, Challenger::new(perm))
}

// #[test]
// fn test_fib() -> Result<(), VerificationError> {
//     type Val = BabyBear;
//     type Challenge = BinomialExtensionField<Val, 4>;

//     type Perm = Poseidon2<Val, DiffusionMatrixBabybear, 16, 7>;
//     let perm = Perm::new_from_rng(8, 22, DiffusionMatrixBabybear, &mut thread_rng());

//     type Dft = Radix2DitParallel;
//     let dft = Dft {};

//     type Challenger = DuplexChallenger<Val, Perm, 16>;

//     let log_n = 0;
//     type Pcs = TrivialPcs<Val, Radix2DitParallel>;
//     let pcs = Pcs::new(log_n, dft);

//     type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
//     let config = MyConfig::new(pcs);
//     let a = 1;
//     let b = 1;
//     let n = 10;
//     let x = 55; 
//     do_test(config, a, b, n, x, Challenger::new(perm))

// }
