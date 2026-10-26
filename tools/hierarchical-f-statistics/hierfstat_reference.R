# Authoring tool: independent reference values from hierfstat (Goudet 2005).
# Usage: Rscript hierfstat_reference.R <counts_csv> <out_json>
# Each locus is analysed on its own with varcomp(); the per-locus
# variance components (region, deme, Ind, Error) are turned into F-statistics.
suppressPackageStartupMessages(library(hierfstat))
args <- commandArgs(trailingOnly = TRUE)
counts <- read.csv(args[1], stringsAsFactors = FALSE)
out <- list()
for (loc in sort(unique(counts$locus))) {
  d <- counts[counts$locus == loc, ]
  region <- c(); deme <- c(); geno <- c()
  for (i in seq_len(nrow(d))) {
    g <- c(rep(11L, d$AA[i]), rep(12L, d$Aa[i]), rep(22L, d$aa[i]))
    region <- c(region, rep(d$region[i], length(g)))
    deme <- c(deme, rep(paste(d$region[i], d$deme[i], sep = "/"), length(g)))
    geno <- c(geno, g)
  }
  vc <- varcomp(data.frame(region = as.integer(factor(region)), deme = as.integer(factor(deme)),
                           L = geno), diploid = TRUE)
  s <- vc$overall               # components summed over the two alleles of the locus
  C <- s[1]; S <- s[2]; I <- s[3]; G <- s[4]
  tot <- C + S + I + G
  out[[loc]] <- list(F_IS = unname(I / (I + G)), F_SC = unname(S / (S + I + G)),
                     F_CT = unname(C / tot), F_IT = unname((C + S + I) / tot),
                     hierfstat_F = vc$F)
  print(loc); print(vc$F)
}
writeLines(jsonlite::toJSON(out, digits = 15, auto_unbox = TRUE, pretty = TRUE), args[2])
