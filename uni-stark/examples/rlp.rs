use std::borrow::Borrow;
use core::mem::{size_of, transmute};
use std::usize;
use itertools::Itertools;
use p3_baby_bear::{BabyBear, DiffusionMatrixBabybear};
use p3_air::{Air, AirBuilder, BaseAir};
use p3_challenger::DuplexChallenger;
use p3_matrix::dense::RowMajorMatrix;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{AbstractExtensionField, AbstractField, ExtensionField, Field, PackedValue, PrimeField, PrimeField32, PrimeField64};
use p3_matrix::{Matrix, MatrixRowSlices};
use p3_poseidon2::Poseidon2;
use p3_uni_stark::{prove, verify, StarkConfig, VerificationError};
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


struct TestPlpCol<T> {
  pub a: [T;3],
  pub b: [T;3],
}

impl<T> Borrow<TestPlpCol<T>> for [T] {
  fn borrow(&self) -> &TestPlpCol<T> {
      debug_assert_eq!(self.len(), NUM_RLP_COLS);
      let (prefix, shorts, suffix) = unsafe { self.align_to::<TestPlpCol<T>>() };
      debug_assert!(prefix.is_empty(), "Alignment should match");
      debug_assert!(suffix.is_empty(), "Alignment should match");
      debug_assert_eq!(shorts.len(), 1);
      &shorts[0]
  }
}

pub struct RlpAir {
  n: u64,
  rlp_array: Vec<u8>,
}

impl Default for RlpAir {
    fn default() -> Self {
        Self {
            n: 0,
            rlp_array: vec![],
        }
    }
}


#[derive(Debug)]
enum RlpValue<'a> {
    String(&'a [u8]),
    List(Vec<RlpValue<'a>>),
}

impl<'a> From<RlpValue<'a>> for Vec<u8> {
  fn from(rlp_value: RlpValue<'a>) -> Self {
      match rlp_value {
          RlpValue::String(bytes) => bytes.to_vec(),
          RlpValue::List(values) => {
              values.into_iter().flat_map(Vec::<u8>::from).collect()
          }
      }
  }
}

const R:u8 = 3;
const RLP_LENGTH:usize = 8;
const RLP_WIDTH:usize = 8;
const MAX_WIDTH:usize = RLP_WIDTH + RLP_LENGTH;

pub(crate) const NUM_RLP_COLS: usize = size_of::<TestPlpCol<u8>>();


impl RlpAir {
    pub fn random_valid_trace<F: PrimeField64>(&self) -> RowMajorMatrix<F>
    {
       
        let rlp_values = decode_rlp_internal(&self.rlp_array).map(|(value, _)| value);
        let mut trace = RowMajorMatrix::new(vec![F::zero(); MAX_WIDTH * 2], MAX_WIDTH);
        
        let rlp_res:Vec::<u64> = match rlp_values {
          Some(RlpValue::String(bytes)) => bytes.into_iter().map(|&x| x as u64).collect(),
          Some(RlpValue::List(bytes)) => bytes.into_iter().map(|rlp_value| {
            match rlp_value {
              RlpValue::String(data) => {
                  let mut res:Vec<u64> = Vec::new();
                  res = data.iter().map(|&x| x as u64).collect();
                  if res.len() < RLP_WIDTH {
                    res.resize(RLP_WIDTH, 0);
                  }

                  let rlc: u64 = data.iter().enumerate().map(|(i, &x)| x as u64 * R.pow(i as u32) as u64).sum();
                  res[RLP_WIDTH-1] = rlc;
                  
                  let mut len = data.len() as u64;
                  let mut binary_len: Vec<u64> = Vec::new();
                  while len > 0 {
                    binary_len.push(len % 2);
                    len /= 2;
                  }
                  if binary_len.len() < RLP_LENGTH {
                    binary_len.resize(RLP_LENGTH, 0);
                  }

                //  let padded_res res.into_iter().chain(binary_len.into_iter());
                let concatenated_array = [res, binary_len].concat(); // 使用 + 运算符连接两个数组
                return concatenated_array;
                  // return res.extend(&binary_len)
              },
              _ => todo!() 
          }
          }).flat_map(Vec::<u64>::from).collect(),
          None => todo!(),
        }; 

        for (i, v) in rlp_res.iter().enumerate() {
          trace.values[i] = F::from_canonical_u64(*v);
        } 

        println!("rlp_res=> {:?}", rlp_res);
        
        trace
    }
}

fn decode_rlp_internal(data: &[u8]) -> Option<(RlpValue, &[u8])> {
match data.get(0)? {
    b if *b <= 0x7f => {
      Some((RlpValue::String(&data[0..1]), &data[1..]))
    }
    b if *b <= 0xb7 => {
      // Short string
      let length = *b as usize - 0x80;
      let value = &data[1..1 + length];
      Some((RlpValue::String(value), &data[1 + length..]))
    }
    b if *b <= 0xbf => {
      // Long string
      let length_bytes = *b as usize - 0xb7;
      let length = decode_length(&data[1..1 + length_bytes])?;
      let value = &data[1 + length_bytes..1 + length_bytes + length];
      Some((RlpValue::String(value), &data[1 + length_bytes + length..]))
    }
    b if *b <= 0xf7 => {
      // Short list
      // let length = *b as usize - 0xc0;
      let mut remaining_data = &data[1..];
      let mut list = Vec::new();
      // for _ in 0..length {
      while !remaining_data.is_empty() {
          let (value, new_remaining_data) = decode_rlp_internal(remaining_data)?;
          list.push(value);
          remaining_data = new_remaining_data;
      }
      Some((RlpValue::List(list), remaining_data))
  }
  b if *b <= 0xff => {
      // Long list
      let length_bytes = *b as usize - 0xf7;
      // let length = decode_length(&data[1..1 + length_bytes])?;
      let mut remaining_data = &data[1 + length_bytes..];
      let mut list = Vec::new();
      while !remaining_data.is_empty() {
          let (value, new_remaining_data) = decode_rlp_internal(remaining_data)?;
          list.push(value);
          remaining_data = new_remaining_data;
      }
      Some((RlpValue::List(list), remaining_data))
  }
    _ => None,
  }
}

fn decode_length(data: &[u8]) -> Option<usize> {
  match data.len() {
      0 => None,
      1 => Some(data[0] as usize),
      _ => {
          let length_bytes = data[0] as usize;
          let mut length = 0;
          for i in 1..=length_bytes {
              length |= (data[i] as usize) << ((length_bytes - i) * 8);
          }
          Some(length)
      }
  }
}


pub fn exp_by_squaring<AF: AbstractField>(val: AF, power_bits: Vec<AF>) -> AF {
  let mut current = val;
  let mut product = AF::one();
  for j in power_bits.into_iter() {
    let c = (current.clone() - AF::one())* j + AF::one();
    product *= c.clone();
    current = current.square();
  }
  product
}

impl<F> BaseAir<F> for RlpAir {
    fn width(&self) -> usize {
      MAX_WIDTH
    }
}

impl<AB: AirBuilder> Air<AB> for RlpAir {
    fn eval(&self, builder: &mut AB) {


        let main = builder.main();
        let local = main.row_slice(0);
        let next = main.row_slice(1);



        let power_bits = [AB::Expr::one(),AB::Expr::zero(),AB::Expr::zero(),AB::Expr::one(),AB::Expr::zero(),AB::Expr::zero(),AB::Expr::zero(),AB::Expr::zero()].to_vec();
        let local_len_arr = local[RLP_WIDTH..].to_vec();
        // let l = local_len_arr.into_iter().map(|x| {
        //   AB::Expr::from_canonical_u8(8);
        // }).collect_vec();

        // AB::Expr::from_canonical_u8(local_len_arr);
        
        let res = exp_by_squaring(AB::Expr::from_canonical_u8(R), power_bits);
        println!("res: {:?}", res);

        let next_len = &next[RLP_WIDTH..];

        // println!("Len: {:?}, {}", next_len, next_len);
        println!("local: {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}", local[0].into(),local[1].into(),local[2].into(),local[3].into(),local[4].into(),local[5].into(),local[6].into(),local[7].into(),local[8].into(),local[9].into(),local[10].into(),local[11].into(),local[12].into(),local[13].into(),local[14].into(),local[15].into());
        println!("next: {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}", next[0].into(),next[1].into(),next[2].into(),next[3].into(),next[4].into(),next[5].into(),next[6].into(),next[7].into(),next[8].into(),next[9].into(),next[10].into(),next[11].into(),next[12].into(),next[13].into(),next[14].into(),next[15].into());
        println!("-------------------------------");

        // let local_acc = local.last();
        // let local_len = local[MAX_WIDTH-2];
        // let next_acc = next.last();
        // let next_len: <AB as AirBuilder>::Var = next[MAX_WIDTH-2];
        
        // let a = AB::F64::as_canonical_u64(next_len);
        // let f = AB::Expr::as_canonical_u64(next_len);
        // assert_eq!(f.as_canonical_u64(), 3);

        // u64::from(next_len);

        // let t =  exp_u64_by_squaring(AB::Expr::from_canonical_u8(R), next_len);
        // let exp_u64 = next_len.into().exp_u64(2);
        // println!("test: {:?}", test);
        // let r = local_len.exp_power_of_2(R as usize);
        // next_acc = next_rlc*R.pow(local_len as u32)
        // let local_len = local[MAX_WIDTH-2];
        // for i in 0..local_len {
        //   local[i]
        // }
        // builder.when_transition()
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
    // 201, 132, 97, 98, 99, 100, 131, 100, 101, 102

    let rlp_array = vec![201, 132, 97, 98, 99, 100, 131, 100, 101, 102];

    println!("rlp_array: {:?}", rlp_array);

    let air = RlpAir { n: 12, rlp_array};

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

    let proof = prove::<MyConfig, _>(&config, &RlpAir::default(), &mut challenger, trace);

    let mut challenger = Challenger::new(perm);
    verify(&config, &RlpAir::default(), &mut challenger, &proof)
}


