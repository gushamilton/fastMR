\name{fast_read_compressed}
\alias{fast_read_compressed}
\title{Read Selected Variants from a CompreSSoR GWAS}
\description{Reads selected columns and canonical variants from a
self-contained Pcodec CompreSSoR store and adds a
\code{chromosome:position:REF:ALT} identity column.}
\usage{
fast_read_compressed(path, variants = NULL,
  columns = c("chromosome", "base_pair_location", "effect_allele",
    "other_allele", "beta", "standard_error"))
}
\arguments{
  \item{path}{Pcodec CompreSSoR store directory.}
  \item{variants}{Optional canonical identity keys or zero-based row IDs.}
  \item{columns}{Columns to return. Identity columns are read internally when
  needed to construct \code{variant_key}.}
}
\value{A data frame containing the requested columns and
\code{variant_key}.}
