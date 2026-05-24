use nalgebra::{DMatrix, DVector};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
  Open,
  Periodic,
}

pub type BasisMatrix = DMatrix<f64>;
pub type SpanMatrix = DMatrix<i64>;
pub type ExtentVector = DVector<i64>;
pub type CoordinateVector = DVector<f64>;
pub type OffsetVector = DVector<i64>;
