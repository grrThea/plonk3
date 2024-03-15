use core::borrow::Borrow;

use alloc::fmt;
use p3_air::{Air, AirBuilder, BaseAir};
use p3_field::AbstractField;
use p3_matrix::MatrixRowSlices;

use crate::columns::{KeccakCols, NUM_KECCAK_COLS};
use crate::constants::rc_value_bit;
use crate::logic::{andn_gen, xor3_gen, xor_gen};
use crate::round_flags::eval_round_flags;
use crate::{BITS_PER_LIMB, NUM_ROUNDS, U64_LIMBS};

pub struct FibAir {}


impl<F> BaseAir<F> for FibAir {
  fn width(&self) -> usize {
      4
  }
}

impl<AB: AirBuilder> Air<AB> for FibAir {
    fn eval(&self, builder: &mut AB) {
      let main = builder.main();
      let _a = main.row_slice(0)[0];
      let _b = main.row_slice(1)[0];
      let _x = main.row_slice(2)[0];
      let _n = main.row_slice(3)[0];

      // println!({}, a);
      // println!({}, b);
      // println!({}, x);
      // println!({}, n);
  
    }
}
