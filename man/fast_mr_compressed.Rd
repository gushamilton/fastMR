\name{fast_mr_compressed}
\alias{fast_mr_compressed}
\title{Run FastMR Directly from Compressed GWAS Files}
\description{Extracts canonical instruments directly from CompreSSoR stores,
reads each outcome once for the union of instruments, and runs every requested
exposure-outcome analysis without rsID lookup or whole-GWAS materialisation.}
\usage{
fast_mr_compressed(exposure_files, outcome_files, instruments,
  methods = "ivw", nboot = 0, seed = NULL, threads = 1,
  io_threads = 1, minimum_snps = 1L, strict = TRUE, output = NULL,
  estimator = c("auto", "pairwise"), steiger = FALSE,
  samplesize_exposure = NULL, samplesize_outcome = NULL,
  steiger_binary = NULL, ...)
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
  \code{minimum_snps}; otherwise omit them with explicit warnings. Each
  warning lists at most \code{getOption("fastMR.warning_pairs", 20)} affected
  pairs/exposures followed by the total count; the per-pair detail is always
  in the \code{counts} element of the \code{compressed_input} attribute. Set
  \code{options(fastMR.warning_pairs = Inf)} for the full listing.}
  \item{output}{Optional path for a Zstandard-compressed Parquet copy of the
  result. The path must not already exist; use \code{\link{fast_write_parquet}}
  when an overwrite or another compression codec is required.}
  \item{estimator}{\code{"auto"} (default) or \code{"pairwise"}. With
  \code{"auto"}, runs requesting only \code{methods = "ivw"} (and no extra
  \code{...} options) whose instrument sets are not shared by all exposures use
  the sparse CSR kernel \code{\link{fast_mr_sparse_ivw}} instead of the
  pairwise path, which assembles one harmonised long table for every pair and
  makes one batched \code{fast_mr()} call. Counts, errors, warnings,
  \code{minimum_snps} handling and row order are identical, but estimates may
  differ from \code{"pairwise"} by rounding (within 1e-14 relative on beta,
  se and Q; usually bit-identical). \code{"pairwise"} and
  every other method set always use the pairwise path. The path used is
  reported as \code{estimator_path} in the \code{compressed_input} attribute.}
  \item{steiger}{If \code{TRUE}, also run per-SNP Steiger filtering
  (\code{\link{fast_mr_steiger_filtering}}) for every retained
  exposure-outcome pair. \code{effect_allele_frequency} and \code{p_value} are
  read in the same pass as beta/standard error, so no store is read twice. The
  MR result itself is unchanged. Steiger rows are the rows MR uses: instruments
  found in both stores with finite beta and positive standard error, for pairs
  meeting \code{minimum_snps}.}
  \item{samplesize_exposure, samplesize_outcome}{Sample sizes for Steiger
  (stores carry none): a vector named by store label (as in
  \code{exposure_files}/\code{outcome_files}; extra names are ignored), an
  unnamed vector in store order, or one value for every store. Required when
  \code{steiger = TRUE}, except for traits listed in \code{steiger_binary}.}
  \item{steiger_binary}{Optional data frame for binary (log-odds) traits with
  columns \code{id} (an exposure or outcome store label), \code{ncase},
  \code{ncontrol} and \code{prevalence}. Those traits use the log-odds
  R-squared model with the harmonic effective sample size, as in
  \code{\link{fast_mr_steiger_filtering}}; a label shared by an exposure and an
  outcome applies to both.}
  \item{...}{Additional options passed to \code{fast_mr()}.}
}
\value{A tidy FastMR result. Extraction counts and input details are retained
in the \code{compressed_input} attribute. With \code{steiger = TRUE}, the
\code{steiger} attribute holds the per-SNP
\code{\link{fast_mr_steiger_filtering}} output (exposure-major, instruments in
supplied order) and \code{compressed_input$timing$steiger_seconds} its compute
time. The \code{steiger} attribute is not written to \code{output}.}
\details{When every exposure has the same instrument set, no requested method
uses bootstrap draws and every pair has all instruments with valid values, the
run uses the shared-grid kernel \code{\link{fast_mr_grid}} whatever
\code{estimator} says (\code{estimator = "pairwise"} does not disable it;
\code{estimator_path} in the \code{compressed_input} attribute is then
\code{"shared_instrument_grid"}). The sparse IVW path is used only when its
estimated memory (native results, the outcome-by-union-instrument matrices it
builds and its per-pair count matrices) fits
\code{getOption("fastMR.sparse_ivw_max_memory_mb", 8192)} MiB; otherwise the
pairwise path runs.

An exposure with an empty instrument set is an error with \code{strict = TRUE}
and is dropped with a warning with \code{strict = FALSE}. A pair with one
matched instrument gets the Wald ratio from the IVW method on every path, as
TwoSampleMR's \code{mr()} reports through \code{mr_wald_ratio()}.}
\seealso{\code{\link{fast_mr}}, \code{\link{fast_read_compressed}},
\code{\link{fast_mr_steiger_filtering}}}
