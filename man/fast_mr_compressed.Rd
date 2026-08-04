\name{fast_mr_compressed}
\alias{fast_mr_compressed}
\title{Run FastMR Directly from Compressed GWAS Files}
\description{Extracts canonical instruments directly from CompreSSoR stores,
reads each outcome once for the union of instruments, and runs every requested
exposure-outcome analysis without rsID lookup or whole-GWAS materialisation.}
\usage{
fast_mr_compressed(exposure_files, outcome_files, instruments,
  methods = "ivw", nboot = 0, seed = NULL, threads = 1,
  io_threads = 1, minimum_snps = 1L, strict = TRUE, ...)
}
\arguments{
  \item{exposure_files}{Named character vector of Pcodec CompreSSoR stores.}
  \item{outcome_files}{Named character vector of Pcodec CompreSSoR stores.}
  \item{instruments}{A canonical-key vector shared by all exposures, or one
  canonical-key vector per exposure in a list.}
  \item{methods}{FastMR method codes.}
  \item{nboot}{Number of bootstrap draws.}
  \item{seed}{Optional FastMR seed.}
  \item{threads}{Native FastMR worker count.}
  \item{io_threads}{Number of stores decoded concurrently inside the shared
  CompreSSoR process.}
  \item{minimum_snps}{Minimum matched instruments required for each pair.}
  \item{strict}{Fail on invalid beta/standard-error values and on pairs below
  \code{minimum_snps}; otherwise omit them with explicit warnings.}
  \item{...}{Additional options passed to \code{fast_mr()}.}
}
\value{A tidy FastMR result. Extraction counts and input details are retained
in the \code{compressed_input} attribute.}
\seealso{\code{\link{fast_mr}}, \code{\link{fast_read_compressed}}}
