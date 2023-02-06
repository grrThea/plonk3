#![no_std]

use hyperfield::field::Field;
use hyperfield::matrix::dense::{DenseMatrix, DenseMatrixView, DenseMatrixViewMut};
use hyperfield::matrix::Matrix;

/// A systematic linear code.
pub trait Code<F: Field> {
    fn encode_batch(&self, messages: DenseMatrixView<F>, codewords: DenseMatrixViewMut<F>);

    fn message_len(&self) -> usize;

    fn codeword_len(&self) -> usize;
}

pub trait SystematicCode<F: Field>: Code<F> {
    fn systematic_len(&self) -> usize;

    fn parity_len(&self) -> usize;

    /// Encode a batch of messages, stored in a matrix with a message in each column.
    ///
    /// Since this is a systemic code, this method extends the input matrix to avoid copying.
    fn append_parity(&self, messages: &mut DenseMatrix<F>) {
        assert_eq!(messages.height(), self.systematic_len());
        messages.expand_to_height(self.codeword_len());
        let (systematic, parity) = messages.as_view_mut().split_rows(self.systematic_len());
        self.write_parity(systematic.as_view(), parity);
    }

    fn write_parity(&self, systematic: DenseMatrixView<F>, parity: DenseMatrixViewMut<F>);
}

impl<F: Field, S: SystematicCode<F>> Code<F> for S {
    fn encode_batch(&self, messages: DenseMatrixView<F>, codewords: DenseMatrixViewMut<F>) {
        let (systematic, parity) = codewords.split_rows(self.systematic_len());
        systematic.values.copy_from_slice(messages.values);
        self.write_parity(messages, parity);
    }

    fn message_len(&self) -> usize {
        self.systematic_len()
    }

    fn codeword_len(&self) -> usize {
        self.systematic_len() + self.parity_len()
    }
}

/// The trivial code whose encoder is the identity function.
pub struct IdentityCode {
    len: usize,
}

impl<F: Field> SystematicCode<F> for IdentityCode {
    fn systematic_len(&self) -> usize {
        self.len
    }

    fn parity_len(&self) -> usize {
        0
    }

    fn write_parity(&self, _systematic: DenseMatrixView<F>, _parity: DenseMatrixViewMut<F>) {
        // All done! There are no parity bits.
    }
}
