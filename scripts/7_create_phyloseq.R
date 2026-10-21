# =============================================================================
#  7. phyloseq objects
#     taxa      : DAS Tool refined MAGs with GTDB-Tk classification (nf-core/mag)
#     functions : KEGG KO, KEGG Module, KEGG Pathway, Pfam, EC, COG (steps 1-6)
#  Abundance = raw estimated reads per sample (not normalised, not rarefied,
#  not transformed):
#     taxa      median MAG depth x MAG length / mean read length
#     functions estreads tables of step 6 (gene depth x gene length / mean read length)
#  Samples   = the samples of the metadata file only (first column = sample ID),
#              all metadata columns (fecal_color, fecal_odor, gender) as sample_data
#  Run:  source("7_create_phyloseq.R")   or   Rscript 7_create_phyloseq.R
# =============================================================================
library(phyloseq)
library(data.table)

## ---- Settings ---------------------------------------------------------------
MAG_OUTDIR      <- "H:/MSU/meg/results"        # nf-core/mag output (refined_bins_only run)
FUNC_DIR        <- "H:/MSU/meg/function"       # functional output (steps 5-6)
SAMPLE_METADATA <- "H:/MSU/meg/metadata.tsv"   # tab-separated; first column = sample ID
OUT_DIR         <- "H:/MSU/meg/phyloseq"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Sample data ------------------------------------------------------------
meta <- fread(SAMPLE_METADATA, colClasses = "character")
samples <- meta[[1]]                                         # sample IDs, in metadata order
if (anyDuplicated(samples)) stop("duplicated sample IDs in metadata: ", paste(unique(samples[duplicated(samples)]), collapse = ", "))
SAM <- sample_data(data.frame(meta, row.names = samples, check.names = FALSE))
in_data <- function(have, what) {
  miss <- setdiff(samples, have)
  if (length(miss)) stop("metadata samples missing from ", what, ": ", paste(miss, collapse = ", "))
}

# mean read length per sample (fastp, after filtering), written by step 5
rl <- fread(file.path(FUNC_DIR, "tables_bacteria", "annotation_summary.tsv"), colClasses = list(character = "sample"))
in_data(rl$sample, "annotation_summary.tsv")
rl <- setNames(rl$mean_read_length, rl$sample)[samples]

## ---- Taxa: GTDB-Tk classified MAGs ------------------------------------------
# bin_depths_summary.tsv : bin, one column per sample = median contig depth of the MAG
# gtdbtk_summary.tsv     : user_genome, classification (d__;p__;c__;o__;f__;g__;s__)
# bin_summary.tsv        : bin, Genome_Size_checkm2 (MAG length, bp)
depth <- fread(file.path(MAG_OUTDIR, "GenomeBinning", "depths", "bins", "bin_depths_summary.tsv"), check.names = FALSE)
gtdb  <- fread(file.path(MAG_OUTDIR, "Taxonomy", "GTDB-Tk", "gtdbtk_summary.tsv"), na.strings = c("", "NA"))
bsum  <- fread(file.path(MAG_OUTDIR, "GenomeBinning", "bin_summary.tsv"), check.names = FALSE)
in_data(names(depth), "bin_depths_summary.tsv")

# MAGs classified by GTDB-Tk that have a depth and a CheckM2 genome size
gtdb <- gtdb[grepl("^d__", classification) & user_genome %in% depth$bin &
             user_genome %in% bsum$bin[!is.na(bsum$Genome_Size_checkm2)]]
if (!nrow(gtdb)) stop("no GTDB-Tk classified MAG matches bin_depths_summary.tsv and bin_summary.tsv")
mags <- gtdb$user_genome

d    <- as.matrix(depth[match(mags, bin), ..samples])
size <- bsum$Genome_Size_checkm2[match(mags, bsum$bin)]
taxa_counts <- sweep(d * size, 2, rl, "/")                   # median depth x MAG size / read length
rownames(taxa_counts) <- sub("\\.fa$", "", mags)
keep <- rowSums(taxa_counts) > 0                             # MAGs absent from all metadata samples
taxa_counts <- taxa_counts[keep, , drop = FALSE]; gtdb <- gtdb[keep]

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# 7 ranks per classification; a rank GTDB-Tk could not assign ("g__", "s__") stays blank ("")
lin <- strsplit(gtdb$classification, ";", fixed = TRUE)
bad <- lengths(lin) != 7
if (any(bad)) stop("classification without 7 ranks: ", paste(head(gtdb$user_genome[bad]), collapse = ", "))
tax <- do.call(rbind, lapply(lin, function(x) trimws(sub("^[dpcofgs]__", "", x))))
dimnames(tax) <- list(rownames(taxa_counts), ranks)

ps_taxa <- phyloseq(otu_table(taxa_counts, taxa_are_rows = TRUE), tax_table(tax), SAM)
stopifnot(!anyNA(otu_table(ps_taxa)), !anyNA(tax_table(ps_taxa)))
saveRDS(ps_taxa, file.path(OUT_DIR, "ps_taxa.rds"))

## ---- Functional groups -----------------------------------------------------
# estreads_<category>.tsv : ID, name (official name, step 6), one column per sample (estimated reads)
categories <- c(KEGG_KO = "KO", KEGG_Module = "KEGG_Module", KEGG_Pathway = "KEGG_Pathway",
                Pfam = "PFAMs", EC = "EC", COG = "COG")
ps_function <- lapply(categories, function(cat) {
  f <- file.path(FUNC_DIR, "tables_bacteria_named", paste0("estreads_", cat, ".tsv"))
  t <- fread(f, colClasses = list(character = 1:2), check.names = FALSE)
  in_data(names(t), basename(f))
  counts <- as.matrix(t[, ..samples]); rownames(counts) <- t[[1]]
  keep <- rowSums(counts) > 0                                         # absent from all metadata samples
  counts <- counts[keep, , drop = FALSE]
  tt <- cbind(ID = t[[1]][keep], NAME = t[[2]][keep]); rownames(tt) <- rownames(counts)   # tax_table: ID, NAME
  ps <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), tax_table(tt), SAM)
  stopifnot(!anyNA(otu_table(ps)), !anyNA(tax_table(ps)))
  ps
})
for (n in names(ps_function)) saveRDS(ps_function[[n]], file.path(OUT_DIR, paste0("ps_", n, ".rds")))
saveRDS(ps_function, file.path(OUT_DIR, "ps_function.rds"))

stopifnot(nsamples(ps_taxa) == length(samples), all(sapply(ps_function, nsamples) == length(samples)))
ps_taxa
ps_function
