\name{fast_mr_masked_ivw}
\alias{fast_mr_masked_ivw}
\title{Low-memory IVW over masked matrices}
\usage{
fast_mr_masked_ivw(
  exposure_beta,
  outcome_beta,
  outcome_se,
  exposure_present,
  outcome_present,
  threads = 1L,
  max_output_cells = 1e+08,
  max_memory_mb = 2048
)
}
\arguments{
\item{exposure_beta}{Numeric \code{E x S} exposure-effect matrix.}
\item{outcome_beta}{Numeric \code{O x S} outcome-effect matrix.}
\item{outcome_se}{Numeric \code{O x S} outcome standard-error matrix.}
\item{exposure_present}{Logical \code{E x S} mask. Values at \code{FALSE}
entries are ignored and may be \code{NA}.}
\item{outcome_present}{Logical \code{O x S} mask. Values at \code{TRUE}
entries must be finite with positive standard errors.}
\item{threads}{Positive integer retained for interface consistency. The dense
path delegates multiplication to the configured BLAS.}
\item{max_output_cells}{Maximum allowed \code{E * O} result cells.}
\item{max_memory_mb}{Conservative bound for native temporary workspace and
result matrices, excluding input matrices already held by R.}
}
\value{
A \code{fastmr_masked_ivw_compact} list containing \code{nsnp}, \code{beta},
\code{se}, \code{Q}, and \code{sigma}, each an \code{E x O} matrix.
}
\description{
Computes all exposure-by-outcome IVW estimates while retaining pair-specific
instrument membership through logical masks. The output and native workspace
are bounded before entering the C++ kernel; process outcomes in batches for
larger scans.
}
\details{
Matrices must already share the same SNP column order. Only entries selected
by the masks are read by the kernel. This function is intended for dense
masked blocks; a sparse CSR exposure panel can avoid dense exposure storage
with \link{fast_mr_sparse_ivw}.
}
