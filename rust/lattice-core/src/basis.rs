use crate::types::BasisMatrix;

#[derive(Clone, Debug, PartialEq)]
pub struct Basis {
  basis: BasisMatrix,
}

impl Basis {
  pub fn new(basis: BasisMatrix) -> Self {
    assert_eq!(basis.nrows(), basis.ncols(), "basis dimension mismatch");
    Self { basis }
  }

  pub fn simple(dim: usize) -> Self {
    Self::new(BasisMatrix::identity(dim, dim))
  }

  pub fn dimension(&self) -> usize {
    self.basis.nrows()
  }

  pub fn basis_vectors(&self) -> &BasisMatrix {
    &self.basis
  }

  pub fn volume(&self) -> f64 {
    self.basis.clone().determinant().abs()
  }
}
