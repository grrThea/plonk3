use itertools::Itertools;
use p3_baby_bear::{BabyBear, DiffusionMatrixBabybear};
use core::mem::{size_of};
use std::arch::x86_64::_CMP_LE_OQ;
use std::clone;

use p3_air::{Air, AirBuilder, BaseAir};
use p3_challenger::DuplexChallenger;
use p3_matrix::dense::RowMajorMatrix;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{Field, PackedValue, PrimeField64};
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

const NUM_RLP_ROWS: usize = 128;

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


impl RlpAir {
    pub fn random_valid_trace<F: PrimeField64>(&self) -> RowMajorMatrix<F>
    {
        let degree = 1;
        let num_rows = NUM_RLP_ROWS * degree;
        let mut trace_values: Vec<F> = vec![F::zero(); num_rows];

        let rlp_values = decode_rlp_internal(&self.rlp_array).map(|(value, _)| value);


        // let a: Vec<_> = match rlp_values {
        //   Some(RlpValue::List(bytes)) => bytes.into_iter().map(|by| {by}).collect_vec(),
        //   Some(RlpValue::String(bytes)) => bytes.to_vec(),
        //   None => todo!(),
        // };

        println!("rlp_values: {:?}", rlp_values);
        let rlp_res = match rlp_values {
          Some(RlpValue::String(bytes)) => bytes.to_vec(),
          Some(RlpValue::List(bytes)) => bytes.into_iter().flat_map(Vec::<u8>::from).collect(),
          None => todo!(),
        }; 

        for (i, v) in rlp_res.iter().enumerate() {

        }

        // let mut j = 0;
        // let remain = vec![0;rlp_res.len()];
        // let values = vec![0;num_rows];
        // while j < 3 * 2 {
        //   let r = remain[j..j+3];
        //   values.into_iter().chain(&remain.into_iter());
        //   j+=3;
        // }

        // let r = vec![0;20];
        // let a = vec![1,2,3,4,5,6,7,8,9];
        // let b = vec![1,2,3,0,4,5,6,0,7,8,9,0,0,0,0,0];

        // let values = vec![0;num_rows];
        // let res = insert_with_padding(values, &rlp_res, 3, 0);
        // println!("res: {:?}", res);

        let mut r = vec![0; 20];
        let a = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let pad_every = 3;
        let pad_value = 0;

        insert_with_padding(&mut r, &a, pad_every, pad_value);

        println!("rs: {:?}", r);
        // for (i, v) in rlp_res.iter().enumerate() {
        //   trace_values[i] = F::from_canonical_u8(*v);
        // } 
        // println!("trace_values: {:?}", trace_values);

        // let mut trace_items: Vec<F> = vec![F::zero();num_rows];
        // for mut i in 0..rlp_res.len(){
        //   trace_items[i] = F::from_canonical_u8(rlp_res[i]);
        //   i+=4;
        // }
        // println!("trace_item: {:?}", trace_items);

        // let trace = vec![vec![F::zero();3];2];
        // let mut trace: Vec<Vec<F>> = vec![vec![F::zero();3];2];
        // for (i, v) in rlp_res.iter().enumerate() {
        //   if i < 3 {
        //     trace[0][i]=F::from_canonical_u8(*v);
        //   } else {
        //     trace[1][i-3]=F::from_canonical_u8(*v);
        //   }
        // } 

        // println!("tracetrace: {:?}", trace);

        RowMajorMatrix::new(trace_values, 4)
    }
}

fn insert_with_padding<T: Copy>(r: &mut Vec<T>, a: &[T], pad_every: usize, pad_value: T) {
  let mut remain = a;
  for (i, &element) in remain.iter().enumerate() {
    r[i] = element;
    if i % pad_every == 0 {
      r[i] = pad_value;
      remain = &a[1+i..];
    }
  }
  // let mut index = 0;
  // // let mut pad_counter = 0;

  // while index < a.len() {
  //     // 插入a中的元素
  //     r[index] = a[index];

  //     // 更新计数器
  //     // pad_counter += 1;
      

  //     // println!("???: {}, {}", index, index % pad_every);
  //     index += 1;
  //     // 检查是否需要插入0
  //     if index % pad_every == 0 && index < a.len() {
  //       println!("indx: {}", index);
  //       r[index] = pad_value;
  //     }
      
  // }
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


impl<F> BaseAir<F> for RlpAir {
    fn width(&self) -> usize {
        4
    }
}

impl<AB: AirBuilder> Air<AB> for RlpAir {
    fn eval(&self, builder: &mut AB) {

        let main: <AB as AirBuilder>::M = builder.main();
        let local = main.row_slice(0);
        let next = main.row_slice(1);
        println!("len: {:?}", local.len());

        println!("local: {:?},{:?},{:?},{:?},{:?},{:?},{:?},{:?}", local[0].into(),local[1].into(),local[2].into(),local[3].into(),local[4].into(),local[5].into(),local[6].into(),local[7].into());
        println!("-------------------------------");
        // builder.when_first_row().assert_one(local[0]);
        // builder.when_first_row().assert_one(local[1]);
        // builder.when_transition().assert_eq(local[0]+local[1], next[0]);
        // builder.when_transition().assert_eq(local[1]+next[0], next[1]);
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

    // "The length of this sentence is more than 55 bytes, I know it because I pre-designed it"
    let rlp_array = vec![200, 131, 97, 98, 99, 131, 100, 101, 102];
    // let rlp_array =vec![248, 94, 131, 97, 98, 99, 248, 88, 179, 84, 104, 101, 32, 108, 101, 110, 103, 116, 104, 32, 111, 102, 32, 116, 104, 105, 115, 32, 115, 101, 110, 116, 101, 110, 99, 101, 32, 105, 115, 32, 109, 111, 114, 101, 32, 116, 104, 97, 110, 32, 53, 53, 32, 98, 121, 116, 101, 115, 44, 32, 163, 73, 32, 107, 110, 111, 119, 32, 105, 116, 32, 98, 101, 99, 97, 117, 115, 101, 32, 73, 32, 112, 114, 101, 45, 100, 101, 115, 105, 103, 110, 101, 100, 32, 105, 116];
    // let rlp_array = vec![248, 88, 179, 84, 104, 101, 32, 108, 101, 110, 103, 116, 104, 32, 111, 102, 32, 116, 104, 105, 115, 32, 115, 101, 110, 116, 101, 110, 99, 101, 32, 105, 115, 32, 109, 111, 114, 101, 32, 116, 104, 97, 110, 32, 53, 53, 32, 98, 121, 116, 101, 115, 44, 32, 163, 73, 32, 107, 110, 111, 119, 32, 105, 116, 32, 98, 101, 99, 97, 117, 115, 101, 32, 73, 32, 112, 114, 101, 45, 100, 101, 115, 105, 103, 110, 101, 100, 32, 105, 116];

    // let rlp_array = vec![184, 86, 84, 104, 101, 32, 108, 101, 110, 103, 116, 104, 32, 111, 102, 32, 116, 104, 105, 115, 32, 115, 101, 110, 116, 101, 110, 99, 101, 32, 105, 115, 32, 109, 111, 114, 101, 32, 116, 104, 97, 110, 32, 53, 53, 32, 98, 121, 116, 101, 115, 44, 32, 73, 32, 107, 110, 111, 119, 32, 105, 116, 32, 98, 101, 99, 97, 117, 115, 101, 32, 73, 32, 112, 114, 101, 45, 100, 101, 115, 105, 103, 110, 101, 100, 32, 105, 116];
    // let rlp_hex: &str = "f902fc0182065c8405f5e1008507081595c083036aab943fc91a3afd70395cd496c647d5a6cc9d4b2b7fad8804fefa17b7240000b902843593564c000000000000000000000000000000000000000000000000000000000000006000000000000000000000000000000000000000000000000000000000000000a000000000000000000000000000000000000000000000000000000000658e953300000000000000000000000000000000000000000000000000000000000000020b080000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002000000000000000000000000000000000000000000000000000000000000004000000000000000000000000000000000000000000000000000000000000000a00000000000000000000000000000000000000000000000000000000000000040000000000000000000000000000000000000000000000000000000000000000200000000000000000000000000000000000000000000000004fefa17b72400000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000004fefa17b72400000000000000000000000000000000000000000000000000251377f5225c886e5e00000000000000000000000000000000000000000000000000000000000000a000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002000000000000000000000000c02aaa39b223fe8d0a0e5c4f27ead9083c756cc2000000000000000000000000269877f972622d3c293fca595c65cf34b7f527cec001a01a8d29640e538156649e5b8d9218a3d91cf866036b907ec5aa31ea69046954e5a002bed5d12542b88d2bb1ae8ba7e838ab75bacddb0aa44703e17de33f0f767ce1";
    // let rlp_array = rlp_hex.to_owned().into_bytes();
    // println!("rlp_array: {:?}", rlp_array);

    // let rlp_array;
    // let test = hex::decode("f902fc0182065c8405f5e1008507081595c083036aab943fc91a3afd70395cd496c647d5a6cc9d4b2b7fad8804fefa17b7240000b902843593564c000000000000000000000000000000000000000000000000000000000000006000000000000000000000000000000000000000000000000000000000000000a000000000000000000000000000000000000000000000000000000000658e953300000000000000000000000000000000000000000000000000000000000000020b080000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002000000000000000000000000000000000000000000000000000000000000004000000000000000000000000000000000000000000000000000000000000000a00000000000000000000000000000000000000000000000000000000000000040000000000000000000000000000000000000000000000000000000000000000200000000000000000000000000000000000000000000000004fefa17b72400000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000004fefa17b72400000000000000000000000000000000000000000000000000251377f5225c886e5e00000000000000000000000000000000000000000000000000000000000000a000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002000000000000000000000000c02aaa39b223fe8d0a0e5c4f27ead9083c756cc2000000000000000000000000269877f972622d3c293fca595c65cf34b7f527cec001a01a8d29640e538156649e5b8d9218a3d91cf866036b907ec5aa31ea69046954e5a002bed5d12542b88d2bb1ae8ba7e838ab75bacddb0aa44703e17de33f0f767ce1");
    // match test {
    //     Ok(data) => rlp_array = data,
    //     Err(message) => panic!("{}", message),
    // };

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


