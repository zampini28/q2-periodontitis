# SMDI computed with the R package microbiome 1.6.0, the implementation
# used by Chen, Marsh & Al-Hebshi (2022), for the cross-check and the
# benchmark run by rcheck.py.
#
#   Rscript smdi_microbiome.R <counts.tsv> <units.tsv> <signature.tsv> <out.tsv>
#   Rscript smdi_microbiome.R --versions
#
# counts.tsv     features x samples read counts (first column: feature id)
# units.tsv      feature id -> unit (species or genus key), as named by
#                q2-periodontitis
# signature.tsv  unit key -> direction (dysbiotic / normobiotic)
#
# Counts are summed per unit, signature units absent from the table are added
# as zero rows, microbiome::transform(x, "clr") is applied to the phyloseq
# object, and SMDI = mean CLR of the dysbiotic units - mean CLR of the
# normobiotic units (NA for an empty sample). Writes one row per sample; the
# last line of stdout is the time in seconds from the loaded tables to the
# index.
suppressMessages({library(microbiome); library(phyloseq)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && args[1] == "--versions") {
  cat(sprintf('{"R": "%s", "microbiome": "%s", "phyloseq": "%s"}\n',
              paste(R.version$major, R.version$minor, sep = "."),
              packageVersion("microbiome"), packageVersion("phyloseq")))
  quit(save = "no")
}
stopifnot(length(args) == 4)
read <- function(path) read.delim(path, check.names = FALSE,
                                  stringsAsFactors = FALSE, quote = "",
                                  comment.char = "")
table <- read(args[1])
counts <- as.matrix(table[, -1, drop = FALSE])
feature_ids <- table[, 1]
units <- read(args[2])
signature <- read(args[3])
unit_of <- units$unit[match(feature_ids, units[, 1])]
stopifnot(!anyNA(unit_of))

seconds <- system.time({
  summed <- rowsum(counts, unit_of, reorder = FALSE)
  absent <- setdiff(signature$unit, rownames(summed))
  if (length(absent)) {
    summed <- rbind(summed, matrix(0, length(absent), ncol(summed),
                                   dimnames = list(absent, colnames(summed))))
  }
  x <- phyloseq(otu_table(summed, taxa_are_rows = TRUE))
  clr <- abundances(microbiome::transform(x, "clr"))
  dys <- signature$unit[signature$direction == "dysbiotic"]
  nor <- signature$unit[signature$direction == "normobiotic"]
  smdi <- colMeans(clr[dys, , drop = FALSE]) -
    colMeans(clr[nor, , drop = FALSE])
  # the CLR of an empty sample is undefined (the plugin gives NaN)
  smdi[colSums(summed) == 0] <- NA
})[["elapsed"]]

write.table(data.frame(sample = names(smdi), SMDI = sprintf("%.17g", smdi)),
            args[4], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%.6f\n", seconds))
