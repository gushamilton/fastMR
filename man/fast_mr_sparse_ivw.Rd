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
  max_memory_mb = 2048,
  pair_snp_keep = NULL
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
\item{pair_snp_keep}{Optional logical \code{O x N} matrix, where \code{N} is
\code{length(col_index)}. Row \code{o} and column \code{k} determine whether
outcome \code{o} keeps the SNP stored at concatenated CSR entry \code{k}.
\code{NA} is rejected; \code{NULL} keeps every entry.}
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
and standard errors at unselected mask entries are ignored. The
\code{pair_snp_keep} columns follow CSR storage order rather than SNP-column
order. For exposure \code{e}, only columns from \code{row_ptr[e] + 1} through
\code{row_ptr[e + 1]} (using R's one-based indexing) are consulted. A mask can
therefore vary by exposure/outcome pair without allocating an
exposure-by-outcome-by-SNP array. The returned \code{nsnp} is the count after
both \code{outcome_present} and \code{pair_snp_keep} are applied.
}
