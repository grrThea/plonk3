use p3_baby_bear::BabyBear;
use p3_challenger::{HashChallenger, SerializingChallenger32};
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_keccak::Keccak256Hash;
use p3_keccak_air::{generate_trace_rows, KeccakAir};
use p3_matrix::Matrix;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_symmetric::{CompressionFunctionFromHasher, SerializingHasher32};
use p3_uni_stark::{prove, verify, StarkConfig, VerificationError};
use p3_util::log2_ceil_usize;
use rand::random;
use tracing_forest::util::LevelFilter;
use tracing_forest::ForestLayer;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};

const NUM_HASHES: usize = 1;
const FIB_WIDTH: usize = 4;

pub struct FibInput {
    pub a: u64,
    pub b: u64,
    pub n: u64,
    pub x: u64,
}

impl FibInput {
    fn new(a: u64, b: u64, x: u64, n: u64) -> Self {
        Self { a, b, x, n }
    }

    fn to_vector(&self) -> Vec<[u64; 25]> {
      vec![[self.a; 25], [self.b; 25],[self.x; 25],[self.n; 25]]
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

    type ByteHash = Keccak256Hash;
    type FieldHash = SerializingHasher32<ByteHash>;
    let byte_hash = ByteHash {};
    let field_hash = FieldHash::new(Keccak256Hash {});

    type MyCompress = CompressionFunctionFromHasher<u8, ByteHash, 2, 32>;
    let compress = MyCompress::new(byte_hash);

    type ValMmcs = FieldMerkleTreeMmcs<Val, u8, FieldHash, MyCompress, 32>;
    let val_mmcs = ValMmcs::new(field_hash, compress);

    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    type Dft = Radix2DitParallel;
    let dft = Dft {};

    type Challenger = SerializingChallenger32<Val, HashChallenger<u8, ByteHash, 32>>;



    let ipts: FibInput = FibInput::new(1, 1, 55, 10);
    let inputs: Vec<[u64; 25]> = ipts.to_vector();
    println!("{:?}", inputs);
    let trace = generate_trace_rows::<Val>(inputs);
    // println!("{:?}", trace);


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

    let mut challenger = Challenger::from_hasher(vec![], byte_hash);

    let proof = prove::<MyConfig, _>(&config, &KeccakAir {}, &mut challenger, trace);

    let mut challenger = Challenger::from_hasher(vec![], byte_hash);
    verify(&config, &KeccakAir {}, &mut challenger, &proof)
}
