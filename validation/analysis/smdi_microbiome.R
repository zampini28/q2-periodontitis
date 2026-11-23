# SMDI computed with the R package microbiome 1.6.0, the implementation
# named by Chen, Marsh & Al-Hebshi (2022, J Dent Res 101:331, Methods:
# "the 'transform' function in R package 'microbiome' (version 1.6.0)"),
# for the cross-check and the benchmark run by rcheck.py.
#
#   Rscript smdi_microbiome.R <counts.tsv> <taxonomy.tsv> <synonyms.tsv> \
#       <signature.tsv> <out.tsv> <units_out.tsv> \
#       [--min-confidence=<c>] [--pseudocount=<p>]
#   Rscript smdi_microbiome.R --units <taxonomy.tsv> <synonyms.tsv> \
#       <level> <units_out.tsv>
#   Rscript smdi_microbiome.R --versions
#
# counts.tsv     features x samples read counts (first column: feature id)
# taxonomy.tsv   feature id -> taxonomy string (column "Taxon"), and
#                optionally the classifier confidence ("Confidence")
# synonyms.tsv   old_name -> canonical_name (q2_periodontitis/data)
# signature.tsv  taxon, direction, level (q2_periodontitis/data)
#
# The taxonomy strings are parsed here, independently of the plugin's code,
# into the unit (species or genus key) each feature is summed into; the
# units are written to units_out.tsv so rcheck.py can compare them with the
# plugin's, feature by feature. Counts are summed per unit, signature units
# absent from the table are added as zero rows, microbiome::transform(x,
# "clr") is applied to the phyloseq object, and SMDI = mean CLR of the
# dysbiotic units - mean CLR of the normobiotic units (NA for an empty
# sample). Writes one row per sample; the last line of stdout is
# "<seconds>\t<peak MB>": the time from the loaded tables and units to the
# index (parsing is not timed, as the plugin's taxonomy parsing is not a
# microbiome function), and the peak RSS of this process (VmHWM, Linux).
#
# Parsing rules (docs/THEORY.md and the docstrings of
# q2_periodontitis/_methods.py):
# * ranks are split on ";". With rank prefixes (g__, s__), the genus and
#   species are the g__ and s__ ranks; without, the 6th and 7th ranks.
#   Placeholder labels (uncultured, unidentified, unclassified, unknown,
#   unassigned, ambiguous, metagenome, environmental, incertae, bacterium,
#   sp) are ignored.
# * names are split on spaces and underscores; eHOMD clades ("clade 431")
#   are dropped, as are database suffix tokens (a number, or one or two
#   capitals: Prevotella_7, Fusobacterium_A) unless they follow "taxon",
#   "hmt", "hot", "ot", "str." or "strain"; a leading "[Genus]" loses its
#   brackets.
# * the genus is the first word, or the first two when the second is an
#   eHOMD placeholder ("Peptostreptococcaceae [G-6]"). A species rank that
#   starts with a capital is a full binomial; otherwise it is an epithet
#   appended to the genus. A species needs at least two words.
# * multi-species labels ("mitis/oralis") are expanded to binomials; such
#   features are not resolved to one species.
# * names are mapped to current names through the synonyms table
#   (case-insensitive).
# * species keys: "hmt N" when the name carries a Human Microbial Taxon
#   number (HMT, HOT or "oral taxon"), the lower-case name otherwise; genus
#   keys: the lower-case genus. Unresolved features keep their reads as
#   "unresolved: <taxonomy string>".
#
# Options (the plugin's --p-min-confidence and --p-smdi-pseudocount):
# * --min-confidence=c: features whose numeric Confidence is below c keep
#   their reads as "low-confidence: <feature id>" and never count as
#   signature taxa; a missing or non-numeric Confidence never excludes.
# * --pseudocount=p: p is added to every relative abundance and the CLR is
#   computed directly, ln(x) - mean(ln(x)) per sample, since microbiome
#   1.6.0 has no fixed-pseudocount option.

read <- function(path) read.delim(path, check.names = FALSE,
                                  stringsAsFactors = FALSE, quote = "",
                                  comment.char = "", na.strings = character(0),
                                  colClasses = "character")

trim <- function(x) gsub("^\\s+|\\s+$", "", x, perl = TRUE)

PLACEHOLDER_STARTS <- c("uncultured", "unidentified", "unclassified",
                        "unknown", "unassigned", "ambiguous", "metagenome",
                        "environmental", "incertae")
ID_TOKENS <- c("taxon", "hmt", "hot", "ot", "str.", "strain")

is_placeholder <- function(name) {
  low <- tolower(trim(name))
  low == "" || any(startsWith(low, PLACEHOLDER_STARTS)) ||
    grepl("metagenome", low, fixed = TRUE) ||
    low %in% c("bacterium", "sp", "sp.")
}

clean_name <- function(raw) {
  raw <- gsub("[\\s_]+clade[\\s_]+\\d+", "", raw, perl = TRUE,
              ignore.case = TRUE)
  tokens <- strsplit(trim(raw), "[\\s_]+", perl = TRUE)[[1]]
  kept <- character(0)
  for (tok in tokens[tokens != ""]) {
    is_suffix <- grepl("^(\\d+|[A-Z]{1,2})$", tok, perl = TRUE)
    if (length(kept) && is_suffix &&
        !(tolower(kept[length(kept)]) %in% ID_TOKENS)) next
    kept <- c(kept, tok)
  }
  sub("^\\[([A-Za-z]+)\\]", "\\1", paste(kept, collapse = " "), perl = TRUE)
}

n_words <- function(x) nchar(gsub("[^ ]", "", x)) + 1

starts_upper <- function(x) grepl("^[[:upper:]]", x, perl = TRUE)

canonical <- function(name, synonyms) {
  hit <- synonyms[tolower(name)]
  unname(ifelse(is.na(hit), name, hit))
}

# "Streptococcus mitis/oralis" -> c("Streptococcus mitis",
# "Streptococcus oralis")
candidates <- function(species) {
  words <- strsplit(species, " ", fixed = TRUE)[[1]]
  parts <- trim(strsplit(paste(words[-1], collapse = " "), "/",
                         fixed = TRUE)[[1]])
  parts <- parts[parts != ""]
  ifelse(starts_upper(parts), parts, paste(words[1], parts))
}

# (genus, species) of one taxonomy string; NA where absent
parse_taxonomy <- function(tax, synonyms) {
  none <- c(genus = NA_character_, species = NA_character_)
  if (is.na(tax) || trim(tax) == "") return(none)
  ranks <- trim(strsplit(tax, ";", fixed = TRUE)[[1]])
  prefixed <- any(grepl("^[a-z]__", ranks, perl = TRUE))
  genus_raw <- NA_character_
  species_raw <- NA_character_
  for (i in seq_along(ranks)) {
    value <- trim(sub("^[a-z]__", "", ranks[i], perl = TRUE))
    if (is_placeholder(value)) next
    if (prefixed) {
      if (startsWith(ranks[i], "g__")) genus_raw <- value
      if (startsWith(ranks[i], "s__")) species_raw <- value
    } else if (i == 6) {
      genus_raw <- value
    } else if (i == 7) {
      species_raw <- value
    }
  }
  genus <- NA_character_
  if (!is.na(genus_raw)) {
    words <- strsplit(clean_name(genus_raw), " ", fixed = TRUE)[[1]]
    if (length(words)) {
      keep <- if (length(words) > 1 && startsWith(words[2], "[")) 2 else 1
      genus <- canonical(paste(words[seq_len(keep)], collapse = " "),
                         synonyms)
    }
  }
  species <- NA_character_
  if (!is.na(species_raw)) {
    epithet <- clean_name(species_raw)
    if (starts_upper(epithet)) {
      species <- epithet
    } else if (!is.na(genus)) {
      species <- paste(genus, epithet)
    }
    if (!is.na(species) && n_words(species) < 2) species <- NA_character_
  }
  if (!is.na(species)) {
    species <- if (grepl("/", species, fixed = TRUE)) {
      paste(canonical(candidates(species), synonyms), collapse = "/")
    } else {
      canonical(species, synonyms)
    }
  }
  c(genus = genus, species = species)
}

HMT <- "(?:\\bHMT|\\bHOT|oral taxon)[\\s_-]*(\\d{1,4})\\b"

taxon_key <- function(name) {
  m <- regmatches(name, regexec(HMT, name, perl = TRUE, ignore.case = TRUE))
  vapply(seq_along(name), function(i) {
    if (length(m[[i]])) sprintf("hmt %d", as.integer(m[[i]][2]))
    else tolower(name[i])
  }, "")
}

unit_of_feature <- function(tax, level, synonyms) {
  name <- parse_taxonomy(tax, synonyms)[[level]]
  if (is.na(name) || name == "" || grepl("/", name, fixed = TRUE)) {
    return(paste0("unresolved: ", if (is.na(tax)) "" else tax))
  }
  if (level == "species") taxon_key(name) else tolower(name)
}

feature_units <- function(feature_ids, taxonomy, level, synonyms,
                          min_confidence = 0) {
  row <- match(feature_ids, taxonomy[, 1])
  taxa <- taxonomy$Taxon[row]
  taxa[is.na(taxa)] <- ""
  units <- vapply(taxa, unit_of_feature, "", level = level,
                  synonyms = synonyms, USE.NAMES = FALSE)
  if ("Confidence" %in% names(taxonomy)) {
    conf <- suppressWarnings(as.numeric(taxonomy$Confidence[row]))
    low <- !is.na(conf) & conf < min_confidence
    units[low] <- paste0("low-confidence: ", feature_ids[low])
  }
  units
}

load_synonyms <- function(path) {
  s <- read(path)
  setNames(trim(s$canonical_name), tolower(trim(s$old_name)))
}

signature_units <- function(signature, synonyms) {
  taxa <- canonical(trim(signature$taxon), synonyms)
  level <- tolower(trim(signature$level[1]))
  list(level = level,
       unit = if (level == "species") taxon_key(taxa) else tolower(taxa),
       direction = tolower(trim(signature$direction)))
}

write_units <- function(feature_ids, units, path) {
  write.table(data.frame(feature = feature_ids, unit = units), path,
              sep = "\t", quote = FALSE, row.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 5 && args[1] == "--units") {
  taxonomy <- read(args[2])
  write_units(taxonomy[, 1],
              feature_units(taxonomy[, 1], taxonomy, args[4],
                            load_synonyms(args[3])), args[5])
  quit(save = "no")
}

suppressMessages({library(microbiome); library(phyloseq)})
if (length(args) == 1 && args[1] == "--versions") {
  cat(sprintf('{"R": "%s", "microbiome": "%s", "phyloseq": "%s"}\n',
              paste(R.version$major, R.version$minor, sep = "."),
              packageVersion("microbiome"), packageVersion("phyloseq")))
  quit(save = "no")
}
option <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit)) as.numeric(sub("^[^=]*=", "", hit[1])) else default
}
min_confidence <- option("min-confidence", 0)
pseudocount <- option("pseudocount", NA)
args <- args[!startsWith(args, "--")]
stopifnot(length(args) == 6)
table <- read(args[1])
counts <- as.matrix(table[, -1, drop = FALSE])
storage.mode(counts) <- "double"
feature_ids <- table[, 1]
taxonomy <- read(args[2])
synonyms <- load_synonyms(args[3])
sig <- signature_units(read(args[4]), synonyms)
units <- feature_units(feature_ids, taxonomy, sig$level, synonyms,
                       min_confidence)

seconds <- system.time({
  summed <- rowsum(counts, units, reorder = FALSE)
  absent <- setdiff(sig$unit, rownames(summed))
  if (length(absent)) {
    summed <- rbind(summed, matrix(0, length(absent), ncol(summed),
                                   dimnames = list(absent, colnames(summed))))
  }
  if (is.na(pseudocount)) {
    x <- phyloseq(otu_table(summed, taxa_are_rows = TRUE))
    clr <- abundances(microbiome::transform(x, "clr"))
  } else {
    logs <- log(sweep(summed, 2, colSums(summed), "/") + pseudocount)
    clr <- sweep(logs, 2, colMeans(logs), "-")
  }
  dys <- sig$unit[sig$direction == "dysbiotic"]
  nor <- sig$unit[sig$direction == "normobiotic"]
  smdi <- colMeans(clr[dys, , drop = FALSE]) -
    colMeans(clr[nor, , drop = FALSE])
  # the CLR of an empty sample is undefined (the plugin gives NaN)
  smdi[colSums(summed) == 0] <- NA
})[["elapsed"]]

write.table(data.frame(sample = names(smdi), SMDI = sprintf("%.17g", smdi)),
            args[5], sep = "\t", quote = FALSE, row.names = FALSE)
write_units(feature_ids, units, args[6])
status <- readLines("/proc/self/status")
peak_kb <- as.numeric(sub("^VmHWM:\\s*([0-9]+) kB$", "\\1",
                          grep("^VmHWM:", status, value = TRUE)))
cat(sprintf("%.6f\t%.3f\n", seconds, peak_kb / 1024))
