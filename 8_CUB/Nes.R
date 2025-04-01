# Functions to calculate the strength of selected codon usage, S = 4*Ne*s
# (2*Ne*s haploids) for amino acids encoded by only two codons.
# Check out the `README' file for some examples.
#
# (c) 2008, Mario dos Reis
#
# If you use this software please cite the appropriate works:
#
# Bulmer (1991) Genetics, 129:897
# [Population genetics model]
#
# Sharp et al. (2005) Nucleic Acids Res. 33:1141
# [First implementation of the model in Eubacteria]
#
# dos Reis and Wernisch (2008) submitted
# [Implementation in Eukaryotes, confidence intervals for NeS]

# A function to calculate S for several genes or group of
# genes simultaneously. Input is usually a codon matrix as the
# one obtained from `codonM'. This matrix might have been
# "compressed" into expresson categories (see README file).
# cf: matrix of codon frequencies, columns are codons
#     and rows are genes/gene sets
# op: optimal codons
# nop: non-optimal codons
# ref: the reference gene set, usually the first row in cf
# mean: whether to return the weighted arithmetic mean for
#       the optimal codons
Nes <- function(cf, op, nop, ref=1, mean=T) {
  cref <- cf[ref,]
  S <- apply(cf, 1, .nesv, flx=cref, op=op, nop=nop, mean=mean)
}

# The most basic function to calculate S for a single two-codon amino acid.
# Phx: relative frequency of optimal codon in highly expressed (HX) genes
# Prf: relative frequency of optimal codon in reference (REF) gene set
.nes <- function(Phx, Prf) {
  return (log(Phx / (1 - Phx)) - log(Prf / (1 - Prf)))
}

# An extension of the previous function to calculate S for several
# amino acids simultaneously.
# fhx: absolute codon frequency vector in HX genes
# flx: absolute codon frequency vector in REF genes
# op: vector indicating which ones are the optimal codons
# nop: vector indicating the corresponding non-optimal codons
# mean: whether to return the weighted average of S for all codon pairs
# alpha: if not NULL, it will calculate the bootstrap percentile confidence
#        interval for the given alpha value
.nesv <- function(fhx, flx, op, nop, mean) {
  # relative frequency of optimal codon:
  Phx <- fhx[op]/(fhx[op] + fhx[nop])
  # relative frequency of nonoptimal codon:
  Prf <- flx[op]/(flx[op] + flx[nop])
  S <- .nes(Phx, Prf)
  if (mean) {
    w <- fhx[op] + fhx[nop]  # weight according to total number of codons
    S <- c(S, Mean=weighted.mean(abs(S), w))
  }
  return(S)
}

# This is a parametric bootstrap approach to calculate the CI for S.
# Fhx: observed number of codons (op and nop) in HX genes
# Frf: observed number of codons (op and nop) in REF genes
# N: number of pseudo bootstrap replicates
# alpha: alpha value for the confidence interval
Nes.boot <- function(Fhx, Frf, N=10000, alpha=0.05) {
  nPhx <- sum(Fhx)
  nPrf <- sum(Frf)
  Phx <- Fhx[1] / nPhx
  Prf <- Frf[1] / nPrf
  
  Phx.boot <- rbinom(N, nPhx, Phx) / nPhx
  Prf.boot <- rbinom(N, nPrf, Prf) / nPrf

  S <- .nes(Phx, Prf)
  S.ci <- .pCI(.nes(Phx.boot, Prf.boot), alpha=alpha)

  return(list(S=S, CI=S.ci))
}

# Percentile confidence interval
.pCI <- function(x, alpha) {
  return(quantile(x, c(alpha/2, 1 - alpha/2), na.rm=T))
}
