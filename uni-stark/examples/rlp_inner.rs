use std::usize;
use p3_baby_bear::{BabyBear, DiffusionMatrixBabybear};
use p3_air::{Air, AirBuilder, BaseAir};
use p3_challenger::DuplexChallenger;
use p3_matrix::dense::RowMajorMatrix;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{AbstractField, Field, PrimeField32, PrimeField64};
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
  rlp_array: Vec<u64>,
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

const PREFIX_NUM:usize = 1;
const LEN_NUM:usize = 3;
const ITEM1_PREFIX_NUM:usize = 1;
const ITEM1_LEN_NUM:usize = 3;
const ITEM1_NUM:usize = 4;
const LOGS_RLP_NUM:usize = 10;
const ITEM2_PREFIX_NUM:usize = 1;
const ITEM2_LEN_NUM:usize = 3;
const ITEM2_NUM:usize = 4;
const FIELD_NUM_ARRAY: [usize; 9] = [PREFIX_NUM, LEN_NUM, ITEM1_PREFIX_NUM, ITEM1_LEN_NUM, ITEM1_NUM, LOGS_RLP_NUM, ITEM2_PREFIX_NUM, ITEM2_LEN_NUM, ITEM2_NUM];

const RLP_FIELD_LEN_NUM:usize = 9;
const RLP_FIELD_NUM: usize = PREFIX_NUM + LEN_NUM + ITEM1_PREFIX_NUM + ITEM1_LEN_NUM + ITEM1_NUM + LOGS_RLP_NUM + ITEM2_PREFIX_NUM + ITEM2_LEN_NUM + ITEM2_NUM;
const TOTAL_FIELD_NUM: usize = RLP_FIELD_NUM + RLP_FIELD_LEN_NUM * 8 + 8 + 1;

// trace field * 20:
// l1_prefix(1), l1_len(3), item1_prefix(1), item1_len(3), item1(4), logs_rlp(10), item2_prefix(1), item2_len(3), item2(4), 
// len(l1_prefix)(1), len(l1_len)(1), len(item1_prefix)(1), len(item1_len)(1), len(item1)(1), len(logs_rlp)(1), len(item2_prefix)(1), len(item2_len)(1), len(item2)(1), total_real_len(1), acc(1)
impl RlpAir {
  pub fn random_valid_trace<F: PrimeField64>(&self) -> RowMajorMatrix<F>
  {

    let data:Vec<u64> = vec![
      207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,343767075299,
      207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,10635747,
      224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,1, 145133,
    ];

    // let data: Vec<u64> = vec![
    //   207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 1, 0, 1, 0, 2, 8, 1, 0, 3, 16, 343767075299,
    //   207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 1, 0, 1, 0, 2, 8, 1, 0, 3, 16, 10635747,
    //   224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 145133,
    // ];
    
    // let data: Vec<u64> = vec![207, 0, 0, 0, 130, 0, 0, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 131, 0, 0, 0, 0, 178, 161, 8, 1, 0, 1, 0, 2, 8, 1, 0, 3, 16, 343767075299,
    // 207, 0, 0, 0, 130, 0, 0, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 131, 0, 0, 0, 0, 178, 161, 8, 1, 0, 1, 0, 2, 8, 1, 0, 3, 16, 10635747,
    // 224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 145133];

    let mut trace = RowMajorMatrix::new(vec![F::zero(); TOTAL_FIELD_NUM * 4], TOTAL_FIELD_NUM);

    for (i, v) in data.iter().enumerate() {
      trace.values[i] = F::from_wrapped_u64(*v);
    } 

    trace
  }     
}


pub fn bits2u64(power_bits: Vec<u64>) -> u64 {
  let num: u64 = power_bits.iter().enumerate().map(|(i, &x)| x as u64 * (2 as u8).pow(i as u32) as u64).sum();
  num
}

pub fn u642bits(len: u8) -> Vec<u8> {
  let mut l = len;
  let mut binary_len: Vec<u8> = Vec::new();
  while l > 0 {
    binary_len.push(len % 2);
    l /= 2;
  }
  binary_len
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

pub fn to_bn<AF: PrimeField64>(val: AF) -> AF {
  val.as_canonical_u64()
}

impl<F> BaseAir<F> for RlpAir {
    fn width(&self) -> usize {
      TOTAL_FIELD_NUM
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
        let r = AB::Expr::from_canonical_u8(R);
        
        let local_acc = local[TOTAL_FIELD_NUM-1];
        let next_acc = next[TOTAL_FIELD_NUM-1];

        // l1_prefix(1), l1_len(3), item1_prefix(1), item1_len(3), item1(4), logs_rlp(10), item2_prefix(1), item2_len(3), item2(4), 
        // len(l1_prefix)(1), len(l1_len)(1), len(item1_prefix)(1), len(item1_len)(1), len(item1)(1), len(logs_rlp)(1), len(item2_prefix)(1), len(item2_len)(1), len(item2)(1), total_real_len(1), acc(1)
        // 207, 0, 0, 0, 130, 0, 0, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 131, 0, 0, 0, 0, 178, 161, 8, 1, 0, 1, 0, 2, 8, 1, 0, 3, 16, 343767075299,
        let mut rlp_arr= vec![AB::Expr::zero(); RLP_FIELD_NUM];
        for i in 0..RLP_FIELD_NUM {
          rlp_arr[i] = local[i].into();
        }
        // println!("rlp_arr: {:?}", rlp_arr);

        let mut rlp_len_arr: Vec<<AB as AirBuilder>::Expr> = vec![];
        for i in RLP_FIELD_NUM..RLP_FIELD_NUM + RLP_FIELD_LEN_NUM * 8 {
          rlp_len_arr.push(local[i].into());
        }
        let nested_vec: Vec<Vec<<AB as AirBuilder>::Expr>> = rlp_len_arr.chunks(8).map(|chunk| chunk.to_vec()).collect();
        // println!("nested_vec: {:?} {:?}", nested_vec.len(), nested_vec);


        // let next_acc_p3: <AB as AirBuilder>::Expr = local_acc + next_exp_sum * exp_u64_by_squaring(AB::Expr::from_canonical_u8(R), local_len_bits);
        let mut acc_p3 = AB::Expr::zero();
        let mut real_len = AB::Expr::zero();
        let mut idx = 0;
        for (i, len) in FIELD_NUM_ARRAY.iter().enumerate() {
          let mut real_len_bits = nested_vec[i].clone();
          let data = rlp_arr[idx..idx + len].to_vec();
      
          let mut rlc = AB::Expr::zero();
          for (i, d) in data.iter().enumerate() {
            rlc += d.clone() * r.exp_u64(i as u64);
          }
          // println!("rlc: {:?}", rlc);
          let acc = acc_p3.clone() + rlc * exp_u64_by_squaring(AB::Expr::from_canonical_u8(R), real_len_bits.clone());
          acc_p3 += acc;
          
          // println!("real_len data: {:?} {:?}", real_len_bits, data);
          idx+=len;
        }

      
        println!("acc_p3: {:?}", acc_p3);
        
        let a = to_bn(AB::Expr::from_canonical_u8(R));
        println!("a=> {:?}", a);

        // let mut local_exp_sum = AB::Expr::zero();
        // for i in 0..TOTAL_FIELD_NUM-1 {
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
    
    
    //rlp: 0xe90282ab34830c781ce0cf82c0afc782c0af83b2a10883b2a108cf82c0afc782c0af83b2a10883b2a108
    //inner: 0xe0cf82c0afc782c0af83b2a10883b2a108cf82c0afc782c0af83b2a10883b2a108
    //inner: [224, 207, 130, 192, 175, 199, 130, 192, 175, 131, 178, 161, 8, 131, 178, 161, 8, 207, 130, 192, 175, 199, 130, 192, 175, 131, 178, 161, 8, 131, 178, 161, 8]
    let hex_string = "e90282ab34830c781ce0cf82c0afc782c0af83b2a10883b2a108cf82c0afc782c0af83b2a10883b2a108";
    let hex_bytes = hex::decode(hex_string).unwrap();

    println!("hex_bytes: {:?}", hex_bytes);

    // let hex_string = "e90282ab34830c781ce0cf82c0afc782c0af83b2a10883b2a108cf82c0afc782c0af83b2a10883b2a108";

    // let rlp_array = hex::decode(hex_string).unwrap();

    // 1,              0,              1,              0,              2,              8,              1,              0,              3
    // 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1
    // 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0
    // let rlp_array:Vec<u64> = vec![
    //   207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 1, 0, 1, 0, 2, 8, 1, 0, 3, 16, 343767075299,
    //   207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 1, 0, 1, 0, 2, 8, 1, 0, 3, 16, 10635747,
    //   224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 145133,
    // ];
    // let rlp_array:Vec<u64> = vec![
    //   207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,343767075299,
    //   207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,10635747,
    //   224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,1, 145133,
    // ];

    let rlp_array:Vec<u64> = vec![
      207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,343767075299,
      207, 0, 0, 0, 130, 0, 0, 0, 192, 175, 0, 0, 199, 130, 192, 175, 131, 178, 161, 8, 0, 0, 131, 0, 0, 0, 178, 161, 8, 0, 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,10635747,
      224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,1, 145133,
    ];


    // let rlp_array = vec![
    //   [224, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 224],
    //   [207, 130, 192, 175, 199, 130, 192, 175, 131, 178, 161, 8, 131, 178, 161, 8, 5245307], 
    //   [207, 130, 192, 175, 199, 130, 192, 175, 131, 178, 161, 8, 131, 178, 161, 8, 5245307],
    // ];
    // println!("rlp_array: {:?}", rlp_array);

    let air = RlpAir { n: 12, rlp_array};

    let trace = air.random_valid_trace::<Val>();
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


