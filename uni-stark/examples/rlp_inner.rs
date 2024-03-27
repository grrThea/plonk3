use std::usize;
use itertools::Itertools;
use p3_baby_bear::{BabyBear, DiffusionMatrixBabybear};
use p3_air::{Air, AirBuilder, BaseAir};
use p3_challenger::DuplexChallenger;
use p3_matrix::dense::RowMajorMatrix;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{AbstractField, Field, PrimeField64};
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

const R:u8 = 2;
const RLP_LENGTH:usize = 8;
const RLP_WIDTH:usize = 8;
const MAX_WIDTH:usize = RLP_WIDTH + RLP_LENGTH;

// outer trace fiels: list prefix(1), list items length(3), item1 
// inner trace fiels: 
impl RlpAir {
  pub fn random_valid_trace<F: PrimeField64>(&self) -> RowMajorMatrix<F>
  {

    // 0xec820acc830bc503e4d183230a88c883230a8883650b3283650b32d182240bc984280a880b83650b3283350e30
    // ["0x0acc", "0x0bc503", [["0x230a88",["0x230a88", "0x650b32"], "0x650b32"], ["0x240b",["0x280a880b", "0x650b32"], "0x350e30"]]]
    // [236, 130, 10, 204, 131, 11, 197, 3, 228, 209, 131, 35, 10, 136, 200, 131, 35, 10, 136, 131, 101, 11, 50, 131, 101, 11, 50, 209, 130, 36, 11, 201, 132, 40, 10, 136, 11, 131, 101, 11, 50, 131, 53, 14, 48]
    // let rlp_values = decode_rlp_internal(&self.rlp_array).map(|(value, _)| value);
    // println!("rlp_values=> {:?}", rlp_values);

    // let rlp_res:Vec::<u64> = match rlp_values {
    //   Some(RlpValue::List(values)) => {
    //     let rlp_arr = values.into_iter().map(|value| {
    //       match value {
    //         // RlpValue::String(data) => {},
    //         RlpValue::List(data) => {
              
    //         },
    //         _ => todo!(),
    //       }
    //     });
    //     vec![]
    //   },
    //   _ => todo!(),
    // };

    
    let trace = RowMajorMatrix::new(vec![F::zero(); MAX_WIDTH * 2], MAX_WIDTH);
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
      // let value = &data[1..1 + length];
      let value = &data[0..1 + length];
      Some((RlpValue::String(value), &data[1 + length..]))
    }
    b if *b <= 0xbf => {
      // Long string
      let length_bytes = *b as usize - 0xb7;
      let length = decode_length(&data[1..1 + length_bytes])?;
      // let value = &data[1 + length_bytes..1 + length_bytes + length];
      let value = &data[0..1 + length_bytes + length];
      Some((RlpValue::String(value), &data[1 + length_bytes + length..]))
    }
    b if *b <= 0xf7 => {
      // Short list
      let mut remaining_data = &data[1..];
      let mut list = Vec::new();
      list.push(RlpValue::String(&data[0..1]));
      while !remaining_data.is_empty() {
          let (value, new_remaining_data) = decode_rlp_internal(remaining_data)?;
          list.push(value);
          remaining_data = new_remaining_data;
      }
      Some((RlpValue::List(list), remaining_data))
  }
  b if *b <= 0xff => {
      // Long list
      let length_bytes: usize = *b as usize - 0xf7;
      // let length = decode_length(&data[1..1 + length_bytes])?;
      let mut remaining_data = &data[1 + length_bytes..];
      let mut list = Vec::new();
      list.push(RlpValue::String(&data[0..1 + length_bytes]));
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

pub fn bits2u64(power_bits: Vec<u64>) -> u64 {
  let num: u64 = power_bits.iter().enumerate().map(|(i, &x)| x as u64 * (2 as u8).pow(i as u32) as u64).sum();
  num
}

pub fn exp_u64_by_squaring<AF: AbstractField>(val: AF, power_bits: Vec<AF>) -> AF {
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

pub fn exp_u64_field_by_squaring<AF: AbstractField>(val: AF, power: u64) -> AF {
  let mut current = val;
  let mut product = AF::one();

  for j in 0..bits_u64(power) {
      if (power >> j & 1) != 0 {
          product *= current.clone();
      }
      current = current.square();
  }
  product
}

const fn bits_u64(n: u64) -> usize {
  (64 - n.leading_zeros()) as usize
}


impl<AB: AirBuilder> Air<AB> for RlpAir {
    fn eval(&self, builder: &mut AB) {

        let main = builder.main();
        let local = main.row_slice(0);
        let next = main.row_slice(1);
        // let r = AB::Expr::from_canonical_u8(R);

        // let local_acc = local[RLP_WIDTH-1];
        // let local_len_bits = local[RLP_WIDTH..].to_vec().iter().map(|&x| x.into()).collect();
        // let next_acc = next[RLP_WIDTH-1];

        // let mut local_exp_sum = AB::Expr::zero();
        // for i in 0..RLP_WIDTH-1 {
        //   local_exp_sum += local[i] * r.exp_u64(i as u64);
        // }

        // let mut next_exp_sum = AB::Expr::zero();
        // for i in 0..RLP_WIDTH-1 {
        //   next_exp_sum += next[i] * r.exp_u64(i as u64);
        // }

        // builder.when_first_row().assert_eq(local_exp_sum, local_acc);

        // let next_acc_p3: <AB as AirBuilder>::Expr = local_acc + next_exp_sum * exp_u64_by_squaring(AB::Expr::from_canonical_u8(R), local_len_bits);
      
        // println!("next_acc: {:?}", next_acc.into());
        // println!("next_acc_p3: {:?}", next_acc_p3);
        // builder.when_transition().assert_eq(next_acc.into(), next_acc_p3);

        // println!("local: {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}", local[0].into(),local[1].into(),local[2].into(),local[3].into(),local[4].into(),local[5].into(),local[6].into(),local[7].into(),local[8].into(),local[9].into(),local[10].into(),local[11].into(),local[12].into(),local[13].into(),local[14].into(),local[15].into());
        // println!("next: {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}", next[0].into(),next[1].into(),next[2].into(),next[3].into(),next[4].into(),next[5].into(),next[6].into(),next[7].into(),next[8].into(),next[9].into(),next[10].into(),next[11].into(),next[12].into(),next[13].into(),next[14].into(),next[15].into());
        // println!("-------------------------------");
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
    
    
    // let hex_string = "ec820acc830bc503e4d183230a88c883230a8883650b3283650b32d182240bc984280a880b83650b3283350e30";

    let hex_string = "e90282ab34830c781ce0cf82c0afc782c0af83b2a10883b2a108cf82c0afc782c0af83b2a10883b2a108";
    // ec820acc830bc503 e4d183230a88c883230a8883650b3283650b32d182240bc984280a880b83650b3283350e30
    // ["0x0acc", "0x0bc503", [["0x230a88",["0x230a88", "0x650b32"], "0x650b32"], ["0x240b",["0x280a880b", "0x650b32"], "0x350e30"]]]
    // [236, 130, 10, 204, 131, 11, 197, 3, 228, 209, 131, 35, 10, 136, 200, 131, 35, 10, 136, 131, 101, 11, 50, 131, 101, 11, 50, 209, 130, 36, 11, 201, 132, 40, 10, 136, 11, 131, 101, 11, 50, 131, 53, 14, 48]
    let rlp_array = hex::decode(hex_string).unwrap();

    println!("rlp_array: {:?}", rlp_array);

    let air = RlpAir { n: 12, rlp_array};

    let trace = air.random_valid_trace();
    println!("trace: {:?}", trace);

    let fri_config = FriConfig {
        log_blowup: 4,
        num_queries: 4,
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


