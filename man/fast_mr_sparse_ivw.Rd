\name{fast_mr_sparse_ivw}
\alias{fast_mr_sparse_ivw}
\title{Low-memory IVW over a sparse CSR exposure panel}
\usage{
fast_mr_sparse_ivw(
  row_ptr,
  col_index,
  exposure_beta,
  outcome_beta,
  outcome_se,
  outcome_present,
  threads = 1L,
  max_output_cells = 1e+08,
  max_memory_mb = 2048
)
}
\arguments{
\item{row_ptr}{Integer CSR row offsets of length \code{E + 1}, starting at
zero and non-decreasing.}
\item{col_index}{Integer zero-based SNP column indices. Indices must be
unique within each exposure row.}
\item{exposure_beta}{Numeric vector of stored exposure effects, parallel to
\code{col_index}.}
\item{outcome_beta}{Numeric \code{O x S} outcome-effect matrix.}
\item{outcome_se}{Numeric \code{O x S} outcome standard-error matrix.}
\item{outcome_present}{Logical \code{O x S} mask. Selected beta values must
be finite with positive standard errors.}
\item{threads}{Positive integer native worker count over exposure rows.}
\item{max_output_cells}{Maximum allowed \code{E * O} result cells.}
\item{max_memory_mb}{Conservative bound for native result workspace,
excluding input matrices already held by R.}
}
\value{
A \code{fastmr_sparse_ivw_compact} list containing \code{nsnp}, \code{beta},
\code{se}, \code{Q}, and \code{sigma}, each an \code{E x O} matrix.
}
\description{
Computes all exposure-by-outcome IVW estimates from a sparse CSR exposure
panel and dense outcome blocks. The output allocation is bounded before the
native kernel is called; process outcomes in batches for larger scans.
}
\details{
CSR column indices are zero-based, unlike ordinary R matrix indices. Outcome
matrices must use the same SNP column order as the CSR columns. Outcome values
and standard errors at unselected mask entries are ignored.
}
