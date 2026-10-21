## =============================================================================
## Shotgun metagenome analysis of the phyloseq objects from 7_create_phyloseq.R
## Adapted from 16S_phyloseq_analysis.R: same steps, functions, statistics,
## filtering rules and figures, applied to every condition of the metadata.
## -----------------------------------------------------------------------------
## Part A - Taxonomy: GTDB-Tk classified MAGs (DAS Tool refined bins, ps_taxa.rds),
##          analysed at every taxonomic level.
## Part B - Function: KEGG KO, KEGG Module, KEGG Pathway, Pfam, EC and COG
##          (ps_function.rds), analysed like the PICRUSt2 part of the 16S
##          workflow; the ggpicrust2 PCA is reproduced with prcomp + ggplot2.
##
## Conditions (metadata columns; reference level first):
##   fecal_color  Normal, Abnormally Dark
##   fecal_odor   Normal, Moderately Unpleasant, Unpleasant
##   gender       female, male
## Every test and condition-coloured figure is produced for each condition
## separately (as for time, Early vs Late, in the 16S workflow). Pairwise tests
## (Wilcoxon, ANCOM-BC2, differential heat tree) are run for every pair of
## levels: one pair for a two-level condition, three for fecal_odor (each on
## the samples of the two levels compared). PERMANOVA / PERMDISP test all
## levels together and, for fecal_odor, every pair as well. Figures of a
## condition show all its levels.
##
## Sections 04-08 (library size, rarefaction, alpha diversity, relative-abundance
## overview, beta diversity) describe the overall community: all bacterial MAGs.
## Sections 09-11 (PERMANOVA, ANCOM-BC2, taxonomic tree) run on each taxonomic
## level, after every level is filtered to exclude unclassified, unknown,
## unnamed and duplicated taxa. PERMANOVA also reports the all-features result.
## Functional tables are filtered in the same way before every analysis.
##
## Main functions: phyloseq, microViz (tax_agg, tax_transform, tax_fix,
## dist_calc, ord_calc, ord_plot, dist_permanova, dist_bdisp, ps_calc_richness,
## ps_calc_diversity, comp_barplot, ps_melt, taxatree_plotkey), microbiome
## (evenness, rarity, prevalence), vegan (rarefy, permutest), ANCOMBC
## (ancombc2), metacoder (heat trees, optional).
##
## Differences from the 16S workflow, required by the data:
##  - one classification (GTDB-Tk): no DADA2 vs QIIME 2 comparison; no Mock sample
##  - no chloroplast / mitochondria filter: GTDB-Tk classifies prokaryotic genomes only
##  - GTDB species names ("Faecalibacterium prausnitzii_C", "CAG-110 sp003524145")
##    are GTDB species clusters; valid when "<genus> <epithet>" with the
##    assigned genus, instead of the SILVA binomial rules
##  - abundances are estimated reads (taxa: median MAG depth x MAG length /
##    mean read length; functions: estimated reads of steps 5-6), rounded to
##    whole reads
##  - rarefaction curves from vegan::rarefy at chosen depths and rarefying by
##    an exact draw without replacement (stats::rhyper; same distribution as
##    phyloseq::rarefy_even_depth(replace = FALSE)): rarecurve(step = 50) and
##    rarefy_even_depth expand every read, impractical for 10^7-10^8 reads
##  - KEGG pathways come from eggNOG-mapper (step 5), so there is no
##    KO -> KEGG conversion; KO features are used for the ANCOM-BC2 bias
##    estimation of KEGG Module and KEGG Pathway (as KO for KEGG in the 16S workflow)
##
## Output layout (one directory per analysis step; condition sub-directories
## where results depend on the condition; ANCOM-BC2 folders per comparison,
## e.g. fecal_odor__Unpleasant_vs_Normal):
##   <out_dir>/00_run_info/
##   <out_dir>/taxa/
##       01_taxonomy  03_filtering  04_library_size  05_alpha_rarefaction
##       06_alpha_diversity  07_relative_abundance  08_beta_diversity
##       09_permanova  10_ancombc2  11_taxonomic_tree
##   <out_dir>/function/
##       02_filtering  03_alpha_diversity  04_beta_diversity  05_permanova  06_ancombc2
##       (01, KO -> KEGG conversion in the 16S workflow, is not needed)
##
## Run:  source("8_metagenome_analysis.R")
##   or  Rscript 8_metagenome_analysis.R --ps_dir H:/MSU/meg/phyloseq --metadata H:/MSU/meg/metadata.tsv
## =============================================================================

## ---- 0. Arguments and configuration -----------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

cfg <- list(
  ps_dir        = get_arg("--ps_dir",   "H:/MSU/meg/phyloseq"),    # ps_taxa.rds, ps_function.rds (7_create_phyloseq.R)
  metadata      = get_arg("--metadata", "H:/MSU/meg/metadata.tsv"),
  mag_dir       = get_arg("--mag_dir",  "H:/MSU/meg/results"),     # nf-core/mag outdir (pipeline_info only)
  out_dir       = get_arg("--out_dir",  "H:/MSU/meg/analysis"),
  threads       = as.integer(get_arg("--threads",
                                     Sys.getenv("SLURM_CPUS_PER_TASK", "1"))),
  seed          = 2026L,
  n_perm        = 9999L,     # permutations for PERMANOVA / PERMDISP
  ancom_prv_cut = 0.10,      # ANCOM-BC2 prevalence filter
  ancom_alpha   = 0.05,
  top_n_abund   = 15L,       # taxa shown individually in relative-abundance plots
  top_n_da      = 5L,        # most significant ANCOM-BC2 features plotted per level
  min_level_reads = 1L,      # samples with no reads at a level are excluded at that level
  rare_points   = 100L,      # depths per rarefaction curve
  rare_saturation = 0.80,    # saturation point: depth at which every sample reaches
                             # this fraction of its observed features (dashed line)
  rare_xmax_factor = 2,      # x axis of rarefaction plots = this x saturation point,
  rare_xmax     = 5e4,       # but at most this many reads (50,000)
  png_dpi       = 300,
  ## Figure styling (all figures: white background)
  font_size     = 15,        # base font size (pt) of every figure
  point_size    = 3.4,       # sample markers
  text_size     = 4.6,       # in-plot labels (mm; 4.6 mm = 13 pt)

  ## Taxon labels treated as unclassified / unknown / unnamed (case-insensitive).
  ## Matched after removal of rank prefixes such as "g__".
  invalid_taxon_regex = paste0(
    "^(na|none|uncultured|unidentified|unknown|unclassified|unassigned|unnamed|",
    "ambiguous|incertae)([^a-z0-9]|$)|metagenome|",
    "(^|[ _])(bacterium|sp\\.?|organism)([ _]|$)"),

  ## Function descriptions treated as unknown / unclassified / unnamed.
  invalid_function_regex = paste0(
    "^(na|none|-|unknown|unclassified|unnamed|uncharacteri[sz]ed|hypothetical)",
    "([^a-z0-9]|$)|",
    "\\b(uncharacteri[sz]ed protein|hypothetical protein|unknown function|",
    "function unknown)\\b"),
  ## 16S rule: features whose description is shared by more than one feature
  ## ID are removed (all copies)
  drop_duplicated_descriptions = TRUE
)

## Conditions: metadata column = its levels, reference level first.
## Two levels: one comparison (second vs first). More levels: every pair of
## levels is compared (later level vs earlier level), each on the samples of
## those two levels; overall tests (PERMANOVA, PERMDISP) use all levels.
conditions <- list(                  # analysed and reported in this order
  gender      = c("female", "male"),
  fecal_color = c("Normal", "Abnormally Dark"),
  fecal_odor  = c("Normal", "Moderately Unpleasant", "Unpleasant"))

## Functional groups: names of the list in ps_function.rds (7_create_phyloseq.R)
fun_groups <- c("KEGG_KO", "KEGG_Module", "KEGG_Pathway", "Pfam", "EC", "COG")

if (is.na(cfg$threads) || cfg$threads < 1) cfg$threads <- 1L

## ---- Packages ----------------------------------------------------------------

required <- c("phyloseq", "microViz", "microbiome", "vegan", "tidyverse", "rstatix",
              "coin", "patchwork", "ANCOMBC")
optional <- c("metacoder", "ggrepel", "yaml", "jsonlite", "writexl", "openxlsx")

missing_req <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_req) > 0) {
  stop("Missing required packages: ", paste(missing_req, collapse = ", "),
       "\nInstall with BiocManager::install(c(",
       paste0('"', missing_req, '"', collapse = ", "), "))",
       "\n(microViz: install.packages(\"microViz\", repos = c(davidbarnett = ",
       "\"https://david-barnett.r-universe.dev\", getOption(\"repos\"))))")
}
has_pkg <- vapply(optional, requireNamespace, logical(1), quietly = TRUE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(microViz)
  library(microbiome)
  library(vegan)
  library(tidyverse)
  library(rstatix)
  library(patchwork)
})

## ---- Helpers -----------------------------------------------------------------

step_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

run_dir  <- step_dir(cfg$out_dir, "00_run_info")
log_file <- file.path(run_dir, "analysis.log")
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
}

## Figure scale factor: fonts, markers, in-plot labels and axis lines are sized
## for a 7 x 7 inch figure and scaled with the figure's size (square root of
## its area relative to 7 x 7), limited to 1-2.5x
fig_scale <- function(width, height) min(max(sqrt(width * height) / 7, 1), 2.5)

## Scales point and text layers of a ggplot by k (markers, p-value labels,
## sample labels); line layers by sqrt(k)
scale_layers <- function(p, k) {
  for (l in p$layers) {
    if (!is.null(l$aes_params$size))      l$aes_params$size      <- l$aes_params$size * k
    if (!is.null(l$aes_params$linewidth)) l$aes_params$linewidth <- l$aes_params$linewidth * sqrt(k)
  }
  p
}

## Writes every figure as PNG on an opaque white background.
##   axes  = TRUE: black x and y axis lines (FALSE for tree layouts)
##   scale = TRUE: font size, markers and labels scaled with the figure size
save_plot <- function(p, dir, name, width = 8, height = 5, dpi = cfg$png_dpi,
                      axes = TRUE, scale = TRUE, k = fig_scale(width, height)) {
  if (!scale) k <- 1
  th <- theme(plot.background = element_rect(fill = "white", colour = NA))
  if (scale) {
    key <- unit(1.25 * 0.9 * cfg$font_size * k, "pt")
    th <- th + theme(text = element_text(size = cfg$font_size * k),
                     legend.key.height = key, legend.key.width = key)
  }
  if (axes) th <- th + theme(axis.line  = element_line(colour = "black", linewidth = 0.7 * sqrt(k)),
                             axis.ticks = element_line(colour = "black", linewidth = 0.6 * sqrt(k)),
                             axis.text  = element_text(colour = "black"))
  if (inherits(p, "patchwork")) {
    if (scale) {
      p <- scale_layers(p, k)
      for (sp in p$patches$plots) if (inherits(sp, "ggplot")) scale_layers(sp, k)
    }
    p <- p & th
  } else if (inherits(p, "ggplot")) {
    if (scale) p <- scale_layers(p, k)
    p <- p + th
  }
  ggsave(file.path(dir, paste0(name, ".png")), p, width = width, height = height,
         dpi = dpi, bg = "white", limitsize = FALSE)
}

save_table <- function(x, dir, name) {
  write_csv(x, file.path(dir, paste0(name, ".csv")))
}

safe_name <- function(x) str_replace_all(tolower(x), "[^a-z0-9]+", "_")

## level text usable in column names: "Abnormally Dark" -> "Abnormally_Dark"
lv_tag <- function(x) str_replace_all(x, "[^A-Za-z0-9]+", "_")

## "fecal_color" -> "Fecal color"
cond_title <- function(cond) str_to_sentence(str_replace_all(cond, "_", " "))

## Write a named list of tables to one Excel workbook, one sheet per table
save_excel <- function(sheets, dir, name) {
  sheets <- sheets[map_lgl(sheets, ~ !is.null(.x) && ncol(.x) > 0)]
  if (length(sheets) == 0) return(invisible(NULL))
  names(sheets) <- str_sub(str_replace_all(names(sheets), "[\\[\\]\\*\\?/\\\\:]", "_"), 1, 31)
  file <- file.path(dir, paste0(name, ".xlsx"))
  if (has_pkg[["writexl"]]) {
    writexl::write_xlsx(map(sheets, as.data.frame), file, format_headers = TRUE)
  } else if (has_pkg[["openxlsx"]]) {
    openxlsx::write.xlsx(map(sheets, as.data.frame), file, overwrite = TRUE,
                         headerStyle = openxlsx::createStyle(textDecoration = "bold"))
  } else {
    log_msg("WARNING: neither writexl nor openxlsx is installed; ", basename(file),
            " not written (the CSV tables contain the same results)")
    return(invisible(NULL))
  }
  log_msg("Excel workbook written: ", file, " (", length(sheets), " sheets)")
  invisible(file)
}

## Outline enclosing all samples of one group: convex hull of the samples with
## rounded corners (hull of circles of radius r around every sample), computed
## in plot-scaled coordinates; r = 3% of the plot span
enclosing_outline <- function(x, y, n = 36, r = 0.03, sx = 1, sy = 1) {
  sx <- if (is.finite(sx) && sx > 0) sx else 1; sy <- if (is.finite(sy) && sy > 0) sy else 1
  xs <- x / sx; ys <- y / sy
  th <- seq(0, 2 * pi, length.out = n + 1)[-1]
  cx <- as.vector(outer(xs, r * cos(th), "+"))
  cy <- as.vector(outer(ys, r * sin(th), "+"))
  h <- grDevices::chull(cx, cy)
  tibble(x = c(cx[h], cx[h[1]]) * sx, y = c(cy[h], cy[h[1]]) * sy)
}

## Silences harmless package notes (microViz short tax_table values,
## ANCOMBC generic matrix input, betadisper negative squared distances for
## the non-Euclidean Bray-Curtis, vegan::rarefy on non-singleton counts:
## estimated reads have few counts of 1; microViz note that the mean relative
## abundances per group are not counts, in the relative-abundance overview)
quiet_pkg_notes <- function(expr) {
  pat <- paste0("Short values detected|tax_fix\\(\\) to make taxa uniquely|",
                "generic 'matrix'/'data.frame' format|some squared distances are negative|",
                "smallest count is|otu_table of counts is NOT available")
  withCallingHandlers(expr,
    message = function(m) if (grepl(pat, conditionMessage(m))) invokeRestart("muffleMessage"),
    warning = function(w) if (grepl(pat, conditionMessage(w))) invokeRestart("muffleWarning"))
}

## Count matrix (samples x features) as a plain numeric matrix
otu_matrix <- function(ps) {
  m <- as(microViz::otu_get(ps), "matrix")
  storage.mode(m) <- "numeric"
  m
}

try_step <- function(what, expr) {
  tryCatch(quiet_pkg_notes(expr), error = function(e) {
    log_msg("WARNING: ", what, " failed: ", conditionMessage(e))
    NULL
  })
}

## Estimated reads rounded to whole reads (rarefying, Chao1 and ANCOM-BC2 need counts)
round_counts <- function(ps) {
  m <- round(as(otu_table(ps), "matrix"))
  otu_table(ps) <- otu_table(m, taxa_are_rows = taxa_are_rows(ps))
  ps
}

cond_cols <- map(conditions, ~ setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")[seq_along(.x)], .x))  # ColorBrewer Dark2
std_ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
alpha_idx <- c("Observed", "Chao1", "Shannon", "Simpson", "Pielou",
               "Low-abundance rarity")

theme_set(theme_bw(base_size = cfg$font_size) +
            theme(panel.grid.minor  = element_blank(),
                  plot.background   = element_rect(fill = "white", colour = NA),
                  panel.background  = element_rect(fill = "white", colour = NA),
                  legend.background = element_rect(fill = "white", colour = NA),
                  legend.key        = element_rect(fill = "white", colour = NA),
                  strip.background  = element_rect(fill = "grey92"),
                  strip.text        = element_text(size = rel(0.95)),
                  plot.title        = element_text(size = rel(1.1), face = "bold"),
                  plot.caption      = element_text(size = rel(0.75), colour = "grey30"),
                  plot.caption.position = "plot",
                  plot.title.position   = "plot",
                  legend.text       = element_text(size = rel(0.9))))

## ---- 1. Inputs ---------------------------------------------------------------

log_msg("Starting analysis; threads = ", cfg$threads)

ps_paths <- c(taxa = file.path(cfg$ps_dir, "ps_taxa.rds"),
              function_groups = file.path(cfg$ps_dir, "ps_function.rds"))
missing_in <- c(ps_paths, metadata = cfg$metadata)
missing_in <- missing_in[!file.exists(missing_in)]
if (length(missing_in) > 0) {
  stop("Input file(s) not found: ", paste(missing_in, collapse = ", "))
}

## Sample design: first column = sample ID; one column per condition.
## The first level of each condition is the reference, so every effect is
## <comparison level> vs <reference level>.
design <- read_tsv(cfg$metadata, col_types = cols(.default = col_character()))
names(design)[1] <- "sample_id"
miss_cond <- setdiff(names(conditions), names(design))
if (length(miss_cond) > 0) {
  stop("Condition column(s) not in ", cfg$metadata, ": ", paste(miss_cond, collapse = ", "),
       "; columns are: ", paste(names(design), collapse = ", "))
}
for (cond in names(conditions)) {
  bad <- setdiff(unique(design[[cond]]), conditions[[cond]])
  if (length(bad) > 0) {
    stop(cond, ": value(s) not in conditions$", cond, " (",
         paste(conditions[[cond]], collapse = ", "), "): ", paste(bad, collapse = ", "))
  }
  design[[cond]] <- factor(design[[cond]], levels = conditions[[cond]])
}
design <- design[str_order(design$sample_id, numeric = TRUE),
                 c("sample_id", names(conditions))]

stopifnot("Duplicated sample IDs in metadata" = !anyDuplicated(design$sample_id))
for (cond in names(conditions)) {
  log_msg("Design, ", cond, ": ",
          paste(names(table(design[[cond]])), table(design[[cond]]), sep = " = ", collapse = ", "))
}
log_msg("Design: ", nrow(design), " samples")
save_table(design, run_dir, "sample_design")

## Pairwise comparisons: one per pair of levels of each condition.
##   cond     metadata column
##   lv       c(reference, comparison); effects are comparison vs reference
##   samples  samples of these two levels
##   id       name used for output folders and files
##   label    "<comparison> vs <reference>"
contrasts <- unlist(lapply(names(conditions), function(cond) {
  lvs <- conditions[[cond]]
  lapply(combn(seq_along(lvs), 2, simplify = FALSE), function(ix) {
    lv <- lvs[ix]
    list(cond = cond, lv = lv,
         samples = design$sample_id[design[[cond]] %in% lv],
         id = if (length(lvs) == 2) cond else paste0(cond, "__", lv_tag(lv[2]), "_vs_", lv_tag(lv[1])),
         label = paste(lv[2], "vs", lv[1]))
  })
}), recursive = FALSE)
names(contrasts) <- map_chr(contrasts, "id")
log_msg("Comparisons: ", paste0(map_chr(contrasts, "cond"), ": ", map_chr(contrasts, "label"),
                                " (n = ", map_int(contrasts, ~ length(.x$samples)), ")", collapse = "; "))
save_table(map_dfr(contrasts, ~ tibble(Comparison_id = .x$id, Condition = .x$cond,
                                       Reference = .x$lv[1], Comparison = .x$lv[2],
                                       n_reference = sum(design[[.x$cond]] == .x$lv[1]),
                                       n_comparison = sum(design[[.x$cond]] == .x$lv[2]))),
           run_dir, "comparisons")

sam_data <- design |>
  as.data.frame() |>
  column_to_rownames("sample_id") |>
  sample_data()

## Restrict a phyloseq object to the design samples and attach the design
attach_design <- function(ps, what) {
  miss <- setdiff(design$sample_id, sample_names(ps))
  if (length(miss) > 0) stop("Design samples missing from ", what, ": ", paste(miss, collapse = ", "))
  ps <- prune_samples(design$sample_id, ps)
  sample_data(ps) <- sam_data
  prune_taxa(taxa_sums(ps) > 0, ps)
}

## Pipeline records (nf-core/mag software versions and applied parameters), if present
record_pipeline_info <- function() {
  info_dir <- file.path(cfg$mag_dir, "pipeline_info")
  out <- list()
  vfile <- list.files(info_dir, pattern = "versions\\.yml$", full.names = TRUE)[1]
  if (!is.na(vfile) && has_pkg[["yaml"]]) {
    sw <- yaml::read_yaml(vfile)
    out$versions <- imap_dfr(sw, function(tools, process) {
      if (!is.list(tools)) return(tibble(Process = process, Tool = process,
                                         Version = as.character(tools)))
      tibble(Process = process, Tool = names(tools),
             Version = map_chr(tools, ~ paste(unlist(.x), collapse = ", ")))
    })
  }
  pfile <- sort(list.files(info_dir, pattern = "^params_.*\\.json$",
                           full.names = TRUE), decreasing = TRUE)[1]
  if (!is.na(pfile) && has_pkg[["jsonlite"]]) {
    p <- jsonlite::read_json(pfile)
    out$params <- tibble(Parameter = names(p),
                         Value = map_chr(p, ~ paste(unlist(.x), collapse = ", ")))
  }
  out
}
pipeline_info <- record_pipeline_info()
if (!is.null(pipeline_info$versions)) save_table(pipeline_info$versions, run_dir, "mag_software_versions")
if (!is.null(pipeline_info$params)) {
  save_table(pipeline_info$params, run_dir, "mag_run_parameters")
  pbi <- pipeline_info$params$Value[pipeline_info$params$Parameter == "postbinning_input"]
  if (length(pbi) == 0 || !pbi[1] %in% "refined_bins_only") {
    log_msg("NOTE: the newest nf-core/mag params file has postbinning_input = ",
            if (length(pbi)) pbi[1] else "<not set>", " (refined_bins_only expected)")
  }
}

## ---- 2. Taxonomy cleaning and level validation -------------------------------

## GTDB-Tk ranks Domain ... Species from 7_create_phyloseq.R (prefixes removed, blank
## where GTDB-Tk assigned no name). Returns the cleaned table plus a per-rank
## log of labels set to missing:
##   missing      - no assignment at that rank
##   placeholder  - label matches invalid_taxon_regex (unclassified/unknown/unnamed)
harmonise_taxonomy <- function(ps_raw) {
  tax_raw <- as.data.frame(as(tax_table(ps_raw), "matrix"), stringsAsFactors = FALSE)
  miss_r <- setdiff(std_ranks, names(tax_raw))
  if (length(miss_r) > 0) stop("Ranks missing from the taxa tax_table: ", paste(miss_r, collapse = ", "))
  tax_mapped <- tax_raw[, std_ranks, drop = FALSE]
  mapping <- tibble(Column = std_ranks, Rank = std_ranks,
                    Source = "GTDB-Tk classification (gtdbtk_summary.tsv)")

  invalid_re <- regex(cfg$invalid_taxon_regex, ignore_case = TRUE)
  stripped   <- mutate(tax_mapped, across(everything(),
                                          ~ str_trim(str_remove(.x, "^\\s*[a-z]__"))))
  raw_labels <- stripped
  rownames(raw_labels) <- rownames(tax_raw)
  is_missing <- map(stripped, ~ is.na(.x) | .x == "")
  is_placeh  <- map(stripped, ~ !is.na(.x) & .x != "" & str_detect(.x, invalid_re))

  placeholders <- imap_dfr(stripped, function(x, r) {
    tibble(Rank = r, Label = x[is_placeh[[r]]])
  }) |>
    count(Rank, Label, name = "Features") |>
    arrange(match(Rank, std_ranks), desc(Features))

  tax_clean <- stripped
  for (r in names(tax_clean)) tax_clean[[r]][is_missing[[r]] | is_placeh[[r]]] <- NA

  ## Known names below an unclassified rank are kept and counted
  inherited <- map(tax_clean, ~ rep(FALSE, length(.x)))
  for (i in seq_along(tax_clean)[-1]) {
    parent_na <- Reduce(`|`, map(tax_clean[seq_len(i - 1)], is.na))
    inherited[[i]] <- parent_na & !is.na(tax_clean[[i]])
  }

  ## Species level: valid GTDB species names only (see valid_species_gtdb())
  vs <- valid_species_gtdb(tax_clean$Genus, tax_clean$Species)
  invalid_sp <- sum(!is.na(tax_clean$Species) & is.na(vs))
  tax_clean$Species <- vs

  label_log <- tibble(
    Rank        = names(tax_clean),
    Features    = nrow(tax_clean),
    Missing     = map_int(names(tax_clean), ~ sum(is_missing[[.x]])),
    Placeholder = map_int(names(tax_clean), ~ sum(is_placeh[[.x]])),
    Known_with_unclassified_parent = map_int(names(tax_clean), ~ sum(inherited[[.x]])),
    Invalid_species_name = if_else(names(tax_clean) == "Species", as.integer(invalid_sp), 0L),
    Valid       = map_int(tax_clean, ~ sum(!is.na(.x))))

  rownames(tax_clean) <- rownames(tax_raw)
  ps <- ps_raw
  tax_table(ps) <- tax_table(as.matrix(tax_clean))

  depth <- tax_clean |>
    summarise(across(everything(), ~ round(100 * mean(!is.na(.x)), 1))) |>
    pivot_longer(everything(), names_to = "Rank", values_to = "Pct_features_classified")

  list(ps = ps, mapping = mapping, depth = depth, label_log = label_log,
       placeholders = placeholders, raw_labels = raw_labels)
}

## Valid GTDB species names: "<genus> <epithet>", two words, the first equal to
## the assigned genus (e.g. "Faecalibacterium prausnitzii_C", "CAG-110 sp003524145")
valid_species_gtdb <- function(genus, species) {
  ok <- !is.na(species) & !is.na(genus) & str_detect(coalesce(species, ""), "^\\S+ \\S+$") &
    word(coalesce(species, ""), 1) == coalesce(genus, "")
  if_else(ok, species, NA_character_)
}

## Duplicated taxa: a name that occurs under conflicting parents at the same
## rank. Processed from the highest to the lowest rank. A feature with an
## unclassified parent is compatible with features of the same name whose
## parent is known, and its missing parent is filled with that consensus value.
## Names with truly conflicting parents are removed at that rank only.
validate_levels <- function(tax) {
  ranks <- intersect(std_ranks, names(tax))
  dup_log <- list(); fill_log <- list()
  for (i in seq_along(ranks)[-1]) {
    r <- ranks[i]; parents <- ranks[1:(i - 1)]
    named <- which(!is.na(tax[[r]]))
    if (length(named) == 0) next
    by_name <- split(named, tax[[r]][named])
    conflict <- names(by_name)[map_lgl(by_name, function(ix) {
      any(map_int(parents, ~ n_distinct(na.omit(tax[ix, .x]))) > 1)
    })]
    if (length(conflict) > 0) {
      lin <- do.call(paste, c(tax[parents], sep = ";"))
      dup_log[[r]] <- tibble(name = tax[[r]], lineage = lin) |>
        dplyr::filter(name %in% conflict) |>
        count(name, lineage, name = "Features") |>
        mutate(Rank = r, .before = 1)
      tax[tax[[r]] %in% conflict & !is.na(tax[[r]]), r] <- NA
    }
    for (nm in setdiff(names(by_name), conflict)) {
      ix <- by_name[[nm]]
      for (pr in parents) {
        v <- unique(na.omit(tax[ix, pr]))
        miss <- ix[is.na(tax[ix, pr])]
        if (length(v) == 1 && length(miss) > 0) {
          tax[miss, pr] <- v
          fill_log[[length(fill_log) + 1]] <- tibble(Rank = r, name = nm, Parent_rank = pr,
                                                     Filled_with = v, Features = length(miss))
        }
      }
    }
  }
  list(tax = tax, duplicates = bind_rows(dup_log), parent_filled = bind_rows(fill_log))
}

## Coverage audit: every label at every rank (Phylum ... Species) found in the
## design samples, with whether its taxon is analysed at that rank and, if not, why
taxa_coverage_audit <- function(raw_labels, clean_tax, design_feats, bac_feats,
                                tax_valid, level_taxa) {
  ranks <- intersect(c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
                     names(raw_labels))
  map_dfr(ranks, function(r) {
    lab <- raw_labels[design_feats, r]
    keep <- !is.na(lab) & str_trim(lab) != ""
    tibble(feature = design_feats[keep], Raw_label = lab[keep]) |>
      group_by(Raw_label) |>
      summarise(features = list(feature), Features = n(), .groups = "drop") |>
      mutate(Rank = r, .before = 1) |>
      mutate(Name = map_chr(features, function(fx) {
               v <- unique(na.omit(clean_tax[fx, r])); if (length(v)) v[1] else NA_character_ }),
             Known = !is.na(Name),
             Analysed = Known & Name %in% level_taxa[[r]],
             Reason = pmap_chr(list(features, Known, Analysed, Name), function(fx, kn, an, nm) {
               if (an) return("")
               if (!kn) return(if (r == "Species") "Not a valid GTDB species name"
                               else "Unclassified / placeholder label")
               if (!any(fx %in% bac_feats)) return("Not Bacteria (Domain)")
               if (all(is.na(tax_valid[intersect(fx, rownames(tax_valid)), r])))
                 return("Same name under conflicting lineages (homonym) - removed")
               if (is.null(level_taxa[[r]])) return("Level not analysed (fewer than 2 taxa or samples)")
               "Not analysed - check"
             }),
             Status = case_when(Analysed ~ "Analysed",
                                Known ~ "Known taxon - excluded",
                                TRUE ~ "Not a known taxon - excluded")) |>
      select(Rank, Raw_label, Name_analysed = Name, Features, Status, Reason)
  }) |>
    mutate(Name_analysed = if_else(Status == "Analysed", Name_analysed, NA_character_))
}

## Level-specific phyloseq object: features with a valid name at `level` are
## kept (phyloseq::prune_taxa) and aggregated per taxon with microViz::tax_agg;
## all other features are excluded. Taxa names = taxon names.
build_level_ps <- function(physeq, tax_valid, level) {
  lvl_ranks <- intersect(std_ranks[1:match(level, std_ranks)], names(tax_valid))
  ps <- physeq
  tax_table(ps) <- tax_table(as.matrix(tax_valid[taxa_names(ps), lvl_ranks, drop = FALSE]))
  keep  <- !is.na(tax_valid[taxa_names(ps), level])
  reads <- taxa_sums(ps)
  ps_k  <- prune_taxa(keep, ps)
  parent_ranks <- setdiff(lvl_ranks, level)
  n_fill <- if (length(parent_ranks)) sum(is.na(as(tax_table(ps_k), "matrix")[, parent_ranks])) else 0L
  if (n_fill > 0) {
    ps_k <- microViz::tax_fix(ps_k, min_length = 1, sep = " ", suffix_rank = "classified",
                              anon_unique = TRUE, verbose = FALSE)
  }
  ps_l  <- microViz::ps_get(microViz::tax_agg(ps_k, rank = level))

  low <- sample_names(ps_l)[sample_sums(ps_l) < cfg$min_level_reads]
  if (length(low) > 0) {
    log_msg("  ", level, ": samples with no reads after filtering excluded at this level: ",
            paste(low, collapse = ", "))
    ps_l <- prune_samples(sample_sums(ps_l) >= cfg$min_level_reads, ps_l)
    ps_l <- prune_taxa(taxa_sums(ps_l) > 0, ps_l)
  }
  list(ps = ps_l,
       log = tibble(Level = level,
                    Features_retained = sum(keep),
                    Features_excluded = sum(!keep),
                    Taxa = ntaxa(ps_l),
                    Parent_labels_filled = n_fill,
                    Reads_retained = sum(reads[keep]),
                    Pct_reads_retained = round(100 * sum(reads[keep]) / sum(reads), 2),
                    Samples = nsamples(ps_l),
                    Samples_excluded = paste(low, collapse = ";")))
}

## ---- 3. Shared analysis helpers ---------------------------------------------

## -- Rarefaction ---------------------------------------------------------------
## Expected number of features (vegan::rarefy) per sample.
## Saturation: for each sample, the smallest depth at which the expected number
## of features reaches cfg$rare_saturation (80%) of its observed number (found
## on a fine log-spaced grid); the saturation point is the largest of these,
## i.e. the depth at which every sample is saturated. The curves are drawn on an
## untransformed read axis from 0 to cfg$rare_xmax_factor x the saturation point,
## at most cfg$rare_xmax reads, at cfg$rare_points evenly spaced depths.
rarefaction_curves <- function(otu_mat) {
  tot <- rowSums(otu_mat)
  rare_at <- function(s, d) as.numeric(quiet_pkg_notes(vegan::rarefy(otu_mat[s, ], sample = d)))
  sat <- map_dfr(rownames(otu_mat), function(s) {
    d <- unique(c(round(10^seq(0, log10(tot[[s]]), length.out = 400)), tot[[s]]))
    f <- rare_at(s, d)
    tibble(sample_id = s, total_reads = tot[[s]], observed = sum(otu_mat[s, ] > 0),
           saturation_depth = d[which(f >= cfg$rare_saturation * sum(otu_mat[s, ] > 0))[1]])
  })
  sat_point <- max(sat$saturation_depth)
  xmax <- min(max(tot), signif(cfg$rare_xmax_factor * sat_point, 2), cfg$rare_xmax, na.rm = TRUE)
  grid <- round(seq(0, xmax, length.out = cfg$rare_points + 1)[-1])
  grid <- sort(unique(c(round(10^seq(0, log10(grid[1]), length.out = 10)), grid)))  # smooth start
  rc <- map_dfr(rownames(otu_mat), function(s) {
    d <- unique(c(grid[grid < tot[[s]]], if (tot[[s]] <= xmax) tot[[s]]))
    tibble(sample_id = s, depth = d, features = rare_at(s, d))
  }) |>
    left_join(design, by = "sample_id")
  attr(rc, "xmax") <- xmax
  attr(rc, "sat_point") <- sat_point
  attr(rc, "saturation") <- sat
  rc
}

plot_rarefaction <- function(rc, cond, rare_depth, title, xlab) {
  xmax <- attr(rc, "xmax"); sat_point <- attr(rc, "sat_point")
  p <- ggplot(rc, aes(depth, features, group = sample_id, colour = .data[[cond]])) +
    geom_line(linewidth = 0.9)
  if (sat_point <= xmax) p <- p +
    geom_vline(xintercept = sat_point, linetype = "dashed", linewidth = 0.7) +
    annotate("label", x = sat_point, y = 0, label = paste0("Saturation: ", scales::comma(sat_point), " reads"),
             hjust = if (sat_point > 0.6 * xmax) 1.05 else -0.05, vjust = 0,
             size = cfg$text_size * 0.75, label.size = 0, fill = "white")
  p +
    scale_colour_manual(values = cond_cols[[cond]]) +
    scale_x_continuous(labels = scales::comma, limits = c(0, xmax), expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = xlab, y = "Expected number of features", colour = cond_title(cond),
         title = title,
         caption = paste0(if (sat_point <= xmax) "Dashed line: s" else "S",
                          "aturation point = ", scales::comma(sat_point),
                          " reads, where every sample reaches ", 100 * cfg$rare_saturation,
                          "% of its observed features",
                          if (sat_point > xmax) " (beyond the plotted range)" else "",
                          ".\nRarefaction depth = ",
                          scales::comma(rare_depth), " reads (smallest sample)",
                          if (rare_depth > xmax) ", beyond the plotted range" else ""))
}

## Rarefying = drawing `depth` reads without replacement from each sample, one
## draw from the multivariate hypergeometric distribution, made feature by
## feature with stats::rhyper (exact; same distribution as
## phyloseq::rarefy_even_depth(replace = FALSE) without expanding every read)
rarefy_counts <- function(m, depth, seed) {
  set.seed(seed)
  out <- t(apply(m, 1, function(x) {
    k <- numeric(length(x)); n_left <- depth; tot_left <- sum(x)
    for (i in which(x > 0)) {
      if (n_left == 0) break
      k[i] <- rhyper(1, x[i], tot_left - x[i], n_left)
      n_left <- n_left - k[i]; tot_left <- tot_left - x[i]
    }
    k
  }))
  dimnames(out) <- dimnames(m)
  out
}

## Rarefied phyloseq object; features left with no reads are removed (trimOTUs)
rarefy_ps <- function(physeq, depth) {
  r <- rarefy_counts(otu_matrix(physeq), depth, cfg$seed)
  r <- r[, colSums(r) > 0, drop = FALSE]
  phyloseq(otu_table(t(r), taxa_are_rows = TRUE), tax_table(physeq), sample_data(physeq))
}

## -- Alpha diversity -----------------------------------------------------------
## Observed and Chao1 richness (microViz::ps_calc_richness), Shannon and
## Gini-Simpson diversity (microViz::ps_calc_diversity), Pielou's evenness
## (microbiome::evenness) and low-abundance rarity (microbiome::rarity)
alpha_values <- function(physeq) {
  ps_a <- physeq |>
    microViz::ps_calc_richness(rank = "unique", index = "observed", varname = "Observed") |>
    microViz::ps_calc_richness(rank = "unique", index = "chao1", varname = "Chao1") |>
    microViz::ps_calc_diversity(rank = "unique", index = "shannon", varname = "Shannon") |>
    microViz::ps_calc_diversity(rank = "unique", index = "gini_simpson", varname = "Simpson")
  even <- microbiome::evenness(physeq, index = "pielou")
  rare <- microbiome::rarity(physeq, index = "low_abundance")
  microViz::samdat_tbl(ps_a, sample_names_col = "sample_id") |>
    select(sample_id, Observed, Chao1, Shannon, Simpson) |>
    mutate(Pielou = even[sample_id, "pielou"],
           `Low-abundance rarity` = rare[sample_id, "low_abundance"]) |>
    mutate(across(all_of(alpha_idx), ~ ifelse(is.finite(.x), .x, NA_real_))) |>
    left_join(design, by = "sample_id") |>
    arrange(match(sample_id, design$sample_id))
}

## Two-sided Wilcoxon rank-sum tests (nominal p-values) for one condition,
## one test per pair of levels (one pair for a two-level condition)
alpha_tests <- function(a, cond) {
  pairs <- keep(contrasts, ~ .x$cond == cond)
  long <- a |>
    mutate(grp = .data[[cond]]) |>
    pivot_longer(all_of(alpha_idx), names_to = "index", values_to = "value") |>
    mutate(index = factor(index, levels = alpha_idx))
  ## indices that do not vary across samples cannot be tested: p = NA
  tests <- long |>
    group_by(index) |>
    group_modify(function(d0, key) map_dfr(pairs, function(ct) {
      lv <- ct$lv
      d <- dplyr::filter(d0, !is.na(value), grp %in% lv) |> mutate(grp = droplevels(grp))
      x <- d$value[d$grp == lv[1]]
      y <- d$value[d$grp == lv[2]]
      base <- tibble(Condition = cond, group1 = lv[1], group2 = lv[2],
                     n1 = length(x), n2 = length(y))
      untestable <- function(note) mutate(base, statistic = NA_real_, p = NA_real_,
                                          effsize = NA_real_, magnitude = NA_character_,
                                          note = note)
      if (length(x) < 2 || length(y) < 2) {
        return(untestable("Fewer than 2 samples with a defined value in a group"))
      }
      if (n_distinct(d$value) < 2) return(untestable("No variation across samples"))
      wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
      es <- tryCatch(rstatix::wilcox_effsize(d, value ~ grp),
                     error = function(e) tibble(effsize = NA_real_, magnitude = NA))
      mutate(base, statistic = unname(wt$statistic), p = wt$p.value,
             effsize = as.numeric(es$effsize), magnitude = as.character(es$magnitude),
             note = "")
    })) |>
    ungroup() |>
    mutate(p.signif = case_when(is.na(p) ~ "NA", p < 0.001 ~ "***", p < 0.01 ~ "**",
                                p < 0.05 ~ "*", TRUE ~ "ns"), .after = p)
  summ <- long |>
    group_by(index, grp) |>
    summarise(Median = median(value, na.rm = TRUE), IQR = IQR(value, na.rm = TRUE),
              n_defined = sum(!is.na(value)), .groups = "drop") |>
    pivot_wider(names_from = grp, values_from = c(Median, IQR, n_defined)) |>
    mutate(Condition = cond, .before = 1)
  list(long = long, tests = tests, summary = summ)
}

## Two levels: the p-value centred above the boxes (as in the 16S workflow).
## More levels: one bracket per pair of levels, labelled with its p-value.
plot_alpha <- function(al, cond, title) {
  rng <- al$long |>
    dplyr::filter(!is.na(value)) |>
    group_by(index) |>
    summarise(top = max(value), span = pmax(diff(range(value)), 1e-9), .groups = "drop")
  ## brackets stacked from the narrowest pair (adjacent levels) upwards
  lab <- left_join(rng, al$tests, by = "index") |>
    mutate(x1 = match(group1, conditions[[cond]]), x2 = match(group2, conditions[[cond]]),
           label = if_else(is.na(p), "not testable", paste0("p = ", signif(p, 2)))) |>
    group_by(index) |>
    arrange(x2 - x1, x1, .by_group = TRUE) |>
    mutate(k = row_number(), y = top + span * (0.18 + 0.30 * (k - 1))) |>
    ungroup()
  n_pairs <- max(lab$k)
  p <- ggplot(dplyr::filter(al$long, !is.na(value)), aes(grp, value, fill = grp)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
    geom_jitter(aes(colour = grp), width = 0.12, height = 0, size = cfg$point_size)
  p <- if (n_pairs == 1) {
    p + geom_text(data = lab, aes(x = 1.5, y = y, label = label),
                  inherit.aes = FALSE, size = cfg$text_size)
  } else {
    p + geom_segment(data = lab, aes(x = x1 + 0.05, xend = x2 - 0.05, y = y, yend = y),
                     inherit.aes = FALSE) +
      geom_text(data = lab, aes(x = (x1 + x2) / 2, y = y, label = label), vjust = -0.35,
                inherit.aes = FALSE, size = cfg$text_size * 0.75)
  }
  p +
    facet_wrap(~ index, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = cond_cols[[cond]]) +
    scale_colour_manual(values = cond_cols[[cond]]) +
    scale_x_discrete(labels = scales::label_wrap(12)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28 + 0.32 * (n_pairs - 1)))) +
    labs(x = NULL, y = "Index value", title = title,
         caption = "p: two-sided Wilcoxon rank-sum test, nominal (unadjusted)") +
    theme(legend.position = "none") +
    ## more than two levels: angled names, so neighbouring labels never overlap
    (if (n_pairs > 1) theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) else NULL)
}

## -- Beta diversity ------------------------------------------------------------
## Bray-Curtis PCoA on relative abundances:
## microViz::tax_transform -> dist_calc -> ord_calc (plotted per condition)
pcoa_bray <- function(physeq) {
  ord <- physeq |>
    microViz::tax_transform("compositional", rank = "unique") |>
    microViz::dist_calc("bray") |>
    microViz::ord_calc(method = "PCoA", verbose = FALSE)
  sc <- as.data.frame(vegan::scores(microViz::ord_get(ord), display = "sites",
                                    choices = 1:2, scaling = 2))
  coords <- sc |> rownames_to_column("sample_id") |> left_join(design, by = "sample_id")
  list(ord = ord, dist = microViz::dist_get(ord), coords = coords, axes = names(sc))
}

## PCoA figure coloured by one condition (microViz::ord_plot), with a rounded
## convex hull around each group and sample labels
plot_pcoa <- function(b, cond, title) {
  ax <- b$axes
  if (length(ax) < 2) {
    log_msg("  PCoA: fewer than 2 ordination axes (too few features) - plot skipped: ", title)
    return(NULL)
  }
  label_layer <- if (has_pkg[["ggrepel"]]) {
    function(...) ggrepel::geom_text_repel(..., seed = cfg$seed, max.overlaps = Inf,
                                           min.segment.length = 0.3, box.padding = 0.4,
                                           point.padding = 0.4)
  } else {
    function(...) geom_text(..., vjust = -1.1)
  }
  p <- microViz::ord_plot(b$ord, axes = 1:2, colour = cond, size = cfg$point_size + 0.6,
                          auto_caption = NA)
  coords <- b$coords
  ## use the plotted coordinates for outlines and labels when they are keyed
  ## by sample name, so that all layers share exactly the same positions
  if (all(ax %in% names(p$data)) && all(rownames(p$data) %in% design$sample_id)) {
    coords <- p$data |>
      as_tibble(rownames = "sample_id") |>
      select(sample_id, all_of(ax)) |>
      left_join(design, by = "sample_id")
  }
  ell <- coords |>
    group_by(.g = .data[[cond]]) |>
    group_modify(~ enclosing_outline(.x[[ax[1]]], .x[[ax[2]]],
                                     sx = diff(range(coords[[ax[1]]])),
                                     sy = diff(range(coords[[ax[2]]])))) |>
    ungroup()
  p +
    geom_polygon(data = ell, aes(x = x, y = y, group = .g, colour = .g, fill = .g),
                 inherit.aes = FALSE, alpha = 0.08, linetype = "dashed", linewidth = 0.9,
                 show.legend = FALSE) +
    label_layer(data = coords, aes(x = .data[[ax[1]]], y = .data[[ax[2]]], label = sample_id,
                                   colour = .data[[cond]]),
                inherit.aes = FALSE, size = cfg$text_size - 0.6, show.legend = FALSE) +
    scale_colour_manual(values = cond_cols[[cond]]) +
    scale_fill_manual(values = cond_cols[[cond]], guide = "none") +
    labs(colour = cond_title(cond), title = title,
         caption = "Dashed outline: all samples of each group (rounded convex hull)")
}

## PCA of a functional profile, reproducing ggpicrust2::pathway_pca (v2.5.19
## source) without the package: features with zero variance removed;
## stats::prcomp(t(abundance), center = TRUE, scale = TRUE); PC1 vs PC2
## points (size 4) coloured by group with ggpicrust2's default colours,
## dashed lines at 0, theme_classic with bold axes; stat_ellipse (default
## t-distribution, level 0.95, white fill, alpha 0.01) for groups with >= 4
## samples; marginal densities of PC1 (top) and PC2 (right), black outline,
## alpha 1, 0.3 of the main plot's height / width. ggpicrust2 joins the panels
## with aplot (axes aligned); here patchwork with the main plot's axis ranges.
plot_pca <- function(physeq, cond, dir, name) {
  abundance <- t(otu_matrix(physeq))                       # features x samples
  abundance <- abundance[apply(abundance, 1, var) != 0, , drop = FALSE]
  pca_result <- stats::prcomp(t(abundance), center = TRUE, scale = TRUE)
  pca_proportion <- (pca_result$sdev[1:2]^2) / sum(pca_result$sdev^2) * 100
  sd <- as(sample_data(physeq), "data.frame")
  pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
                         Group = droplevels(factor(sd[rownames(pca_result$x), cond])),
                         row.names = rownames(pca_result$x))
  colors <- setNames(c("#d93c3e", "#3685bc", "#208A42", "#89288F")[seq_len(nlevels(pca_data$Group))],
                     levels(pca_data$Group))                # ggpicrust2 default colours
  pca_plot <- ggplot(pca_data, aes(x = .data$PC1, y = .data$PC2)) +
    geom_point(size = 4, aes(color = .data$Group), show.legend = TRUE) +
    scale_color_manual(values = colors) +
    labs(x = paste0("PC1(", round(pca_proportion[1], 1), "%)"),
         y = paste0("PC2(", round(pca_proportion[2], 1), "%)"), color = cond) +
    theme_classic() +
    theme(axis.line = element_line(colour = "black"),
          axis.title = element_text(color = "black", face = "bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(color = "black", size = 10, face = "bold"),
          legend.text = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 16, face = "bold")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  ellipse_groups <- names(which(table(pca_data$Group) >= 4))
  if (length(ellipse_groups) > 0) {
    pca_plot <- pca_plot +
      stat_ellipse(data = pca_data[pca_data$Group %in% ellipse_groups, , drop = FALSE],
                   aes(color = .data$Group), fill = "white", geom = "polygon",
                   level = 0.95, alpha = 0.01, show.legend = FALSE)
  }
  ## marginal densities aligned to the main plot's axis ranges
  pp <- ggplot_build(pca_plot)$layout$panel_params[[1]]
  xr <- pp$x.range; yr <- pp$y.range
  pca_plot <- pca_plot + coord_cartesian(xlim = xr, ylim = yr, expand = FALSE)
  density_data <- pca_data[pca_data$Group %in% names(which(table(pca_data$Group) >= 2)), , drop = FALSE]
  pc1_density <- ggplot(density_data) +
    geom_density(aes(x = .data$PC1, group = .data$Group, fill = .data$Group),
                 color = "black", alpha = 1, position = "identity", show.legend = FALSE) +
    scale_fill_manual(values = colors) + theme_classic() +
    scale_y_continuous(expand = c(0, 0.001)) + labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    coord_cartesian(xlim = xr, expand = FALSE)
  pc2_density <- ggplot(density_data) +
    geom_density(aes(x = .data$PC2, group = .data$Group, fill = .data$Group),
                 color = "black", alpha = 1, position = "identity", show.legend = FALSE) +
    scale_fill_manual(values = colors) + theme_classic() +
    scale_y_continuous(expand = c(0, 0.001)) + labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    coord_flip(xlim = yr, expand = FALSE)
  p <- pc1_density + plot_spacer() + pca_plot + pc2_density +
    plot_layout(ncol = 2, widths = c(1, 0.3), heights = c(0.3, 1), guides = "collect")
  ## the ggpicrust2 figure keeps its own sizes (in the 16S workflow it is a
  ## ggplotify grob, so save_plot's scaling did not reach it)
  save_plot(p, dir, name, 9, 8, axes = FALSE, scale = FALSE)
}

## PERMANOVA (microViz::dist_permanova, vegan::adonis2) and homogeneity of
## dispersions (microViz::dist_bdisp, vegan::betadisper; permutation test with
## vegan::permutest) on Bray-Curtis distances of relative abundances, one condition.
## ct = NULL: all levels of the condition ("All groups"); ct = a pairwise
## comparison: the samples of its two levels only.
permanova_bray <- function(physeq, label, cond, raw_dir = NULL, ct = NULL) {
  comparison <- "All groups"
  if (!is.null(ct)) {
    physeq <- prune_samples(intersect(ct$samples, sample_names(physeq)), physeq)
    physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
    comparison <- ct$label
  }
  sd <- as(sample_data(physeq), "data.frame")
  sd[[cond]] <- droplevels(factor(sd[[cond]], levels = conditions[[cond]]))
  sample_data(physeq) <- sample_data(sd)
  d <- physeq |>
    microViz::tax_transform("compositional", rank = "unique") |>
    microViz::dist_calc("bray")
  perm <- microViz::perm_get(
    microViz::dist_permanova(d, variables = cond, n_perms = cfg$n_perm,
                             seed = cfg$seed, by = "terms", verbose = FALSE))
  bd <- microViz::bdisp_get(
    microViz::dist_bdisp(d, variables = cond, method = "centroid", verbose = FALSE))
  set.seed(cfg$seed)
  pd <- vegan::permutest(bd[[cond]]$model, permutations = cfg$n_perm)$tab

  ## Raw outputs, one set of TSV files per level / functional group:
  ##   permanova_raw_<x>.tsv   - complete adonis2 table (condition, Residual, Total)
  ##   permdisp_raw_<x>.tsv    - complete permutest table of the betadisper model
  ##   permdisp_distances_<x>.tsv - distance of each sample to its group centroid
  if (!is.null(raw_dir)) {
    nm <- safe_name(if (is.null(ct)) label else paste(label, ct$label))
    hdr <- function(x) mutate(x, Condition = cond, Comparison = comparison, Level = label,
                              Distance = "Bray-Curtis",
                              Permutations = cfg$n_perm, Seed = cfg$seed, .before = 1)
    write_tsv(hdr(as_tibble(as.data.frame(perm), rownames = "Term")),
              file.path(raw_dir, paste0("permanova_raw_", nm, ".tsv")), na = "")
    write_tsv(hdr(as_tibble(as.data.frame(pd), rownames = "Term")),
              file.path(raw_dir, paste0("permdisp_raw_", nm, ".tsv")), na = "")
    write_tsv(tibble(sample_id = names(bd[[cond]]$model$distances),
                     group = as.character(bd[[cond]]$model$group),
                     distance_to_centroid = unname(bd[[cond]]$model$distances)),
              file.path(raw_dir, paste0("permdisp_distances_", nm, ".tsv")))
  }
  tibble(Condition = cond, Comparison = comparison, Level = label, Features = ntaxa(physeq),
         F = perm[cond, "F"], R2 = perm[cond, "R2"], p = perm[cond, "Pr(>F)"],
         PERMDISP_F = pd$F[1], PERMDISP_p = pd$`Pr(>F)`[1])
}

## ---- ANCOM-BC2 ------------------------------------------------------------
## ANCOM-BC2 (Lin & Peddada 2024), <comparison level> vs <reference level> of
## one condition, on every taxon of a level / feature of a functional group.
##
## Taxonomy check before testing (verify_level_taxonomy): every rank column of
## the level's tax_table (Domain ... level) must carry a name. Missing,
## empty, placeholder or ambiguous labels are re-labelled with
## microViz::tax_fix(), and every change is logged.
##
## Two ANCOM-BC2 fits per table, so that every feature receives a result:
##  1. primary fit - zeros treated as missing on the log scale (ANCOM-BC2
##     default, no pseudo-count), prevalence filter cfg$ancom_prv_cut.
##     For the taxonomic levels the sample-specific biases are estimated from
##     all bacterial MAGs (`data`) and the tests run on the level's taxa
##     (`aggregate_data`); for KEGG Module / Pathway the biases are estimated
##     from the KO features. Features whose variance is zero are removed from
##     the bias estimation only (named in the package's "Zero variances" error).
##  2. pseudo-count fit - for features the primary fit cannot estimate: a
##     pseudo-count of 1 is added to all counts and no prevalence filter is
##     applied, so every feature can be estimated.
## q-values are Holm-adjusted over all features of the table. Every feature
## also gets a two-sided Wilcoxon rank-sum test of relative abundance (nominal p).
## Features are passed under neutral IDs (A00001 for bias features, F00001 for
## tested features); the condition is passed as `group` (Ref / Test), so that
## names with spaces or brackets cannot be altered.

## Verify, and re-label where necessary, the taxon names of a level object
verify_level_taxonomy <- function(physeq, level) {
  ranks <- rank_names(physeq)
  tt0 <- as.data.frame(as(tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
  invalid_re <- regex(cfg$invalid_taxon_regex, ignore_case = TRUE)
  check <- function(tt) {
    map_dfr(seq_along(ranks), function(i) {
      r <- ranks[i]; x <- tt[[r]]
      parent <- if (i > 1) do.call(paste, c(tt[ranks[1:(i - 1)]], sep = ";")) else ""
      ok <- !is.na(x) & str_trim(x) != ""
      dup <- tibble(x = x, parent = parent) |> dplyr::filter(!is.na(x)) |>
        distinct() |> count(x) |> dplyr::filter(n > 1) |> pull(x)
      tibble(Level = level, Rank = r, Taxa = nrow(tt),
             Missing = sum(!ok),
             Placeholder = sum(ok & str_detect(x, invalid_re)),
             Ambiguous = sum(x %in% dup))
    })
  }
  need_ranks <- std_ranks[seq_len(match(level, std_ranks))]
  for (r in setdiff(need_ranks, ranks)) tt0[[r]] <- NA_character_
  ranks <- c(need_ranks, setdiff(ranks, need_ranks))
  tt0 <- tt0[, ranks, drop = FALSE]
  if (!identical(rank_names(physeq), ranks)) tax_table(physeq) <- tax_table(as.matrix(tt0))
  before <- check(tt0)
  needs_fix <- any(before$Missing > 0 | before$Placeholder > 0 | before$Ambiguous > 0)
  relabel <- tibble(Level = character(), feature = character(), Rank = character(),
                    Original = character(), Relabelled = character())
  if (needs_fix) {
    tt_na <- tt0
    for (r in ranks) {
      bad <- is.na(tt_na[[r]]) | str_trim(tt_na[[r]]) == "" |
        str_detect(coalesce(tt_na[[r]], ""), invalid_re)
      tt_na[[r]][bad] <- NA
    }
    tax_table(physeq) <- tax_table(as.matrix(tt_na))
    physeq <- microViz::tax_fix(physeq, min_length = 1, sep = " ",
                                suffix_rank = "classified", anon_unique = TRUE,
                                verbose = FALSE)
    tt1 <- as.data.frame(as(tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
    relabel <- map_dfr(ranks, function(r) {
      ch <- coalesce(tt0[[r]], "<missing>") != tt1[[r]]
      tibble(Level = level, feature = rownames(tt0)[ch], Rank = r,
             Original = tt0[[r]][ch], Relabelled = tt1[[r]][ch])
    })
    taxa_names(physeq) <- make.unique(tt1[[level]], sep = " ")
  }
  after <- check(as.data.frame(as(tax_table(physeq), "matrix"), stringsAsFactors = FALSE))
  report <- left_join(before, select(after, Rank, Missing_after = Missing,
                                     Placeholder_after = Placeholder,
                                     Ambiguous_after = Ambiguous), by = "Rank") |>
    mutate(Relabelled = map_int(Rank, ~ sum(relabel$Rank == .x)))
  list(ps = physeq, report = report, relabel = relabel)
}

## Per-feature input for one condition: count matrix with neutral IDs,
## metadata, descriptive statistics and Wilcoxon tests
##   prevalence        microbiome::prevalence (samples with count > 0)
##   relative abund.   microViz::tax_transform("compositional")
prepare_da_input <- function(physeq, ct) {
  cond <- ct$cond; lv <- ct$lv; tg <- lv_tag(lv)
  physeq <- prune_samples(intersect(ct$samples, sample_names(physeq)), physeq)   # the two levels
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  mat <- t(otu_matrix(physeq))                          # features x samples (counts)
  ids <- sprintf("F%05d", seq_len(nrow(mat)))
  id_map <- setNames(rownames(mat), ids)
  rownames(mat) <- ids
  sd <- as(sample_data(physeq), "data.frame")[colnames(mat), , drop = FALSE]
  ## sample_id is kept as a second column: ANCOMBC's data check subsets
  ## meta_data[samples, ], which drops a one-column data frame to a vector
  meta_data <- data.frame(sample_id = colnames(mat),
                          group = factor(if_else(sd[[cond]] == lv[2], "Test", "Ref"),
                                         levels = c("Ref", "Test")),
                          row.names = colnames(mat))

  ps_rel <- microViz::ps_get(microViz::tax_transform(physeq, "compositional", rank = "unique"))
  s_by   <- map(set_names(lv), ~ sample_names(physeq)[sd[sample_names(physeq), cond] == .x])
  ps_by  <- map(s_by, ~ prune_samples(.x, physeq))
  rel_by <- map(s_by, ~ prune_samples(.x, ps_rel))
  feat    <- unname(id_map)
  prev_n  <- map(ps_by, ~ microbiome::prevalence(.x, detection = 0, count = TRUE)[feat])
  mean_ra <- map(rel_by, ~ rowMeans(microbiome::abundances(.x))[feat] * 100)
  rel_all <- microbiome::abundances(ps_rel)[feat, , drop = FALSE] * 100
  is1 <- colnames(rel_all) %in% s_by[[1]]
  is2 <- colnames(rel_all) %in% s_by[[2]]
  wilcox_p <- apply(rel_all, 1, function(x) {
    if (n_distinct(x) < 2) return(NA_real_)
    suppressWarnings(wilcox.test(x[is1], x[is2], exact = FALSE)$p.value)
  })
  const_in <- function(p) apply(microbiome::abundances(p)[feat, , drop = FALSE], 1,
                                function(x) n_distinct(x) <= 1)
  c1 <- const_in(ps_by[[1]]); c2 <- const_in(ps_by[[2]])
  desc_stats <- tibble(
    id    = ids,
    taxon = feat,
    prevalence_all = round((prev_n[[1]] + prev_n[[2]]) / nsamples(physeq), 3),
    !!paste0("n_present_", tg[1]) := unname(prev_n[[1]]),
    !!paste0("n_present_", tg[2]) := unname(prev_n[[2]]),
    !!paste0("mean_rel_abund_", tg[1]) := unname(mean_ra[[1]]),
    !!paste0("mean_rel_abund_", tg[2]) := unname(mean_ra[[2]]),
    wilcox_p_nominal = unname(wilcox_p),
    constant_in = case_when(c1 & c2 ~ paste(lv[1], "and", lv[2]), c1 ~ lv[1],
                            c2 ~ lv[2], TRUE ~ NA_character_),
    absent_in = case_when(prev_n[[1]] == 0 ~ lv[1],
                          prev_n[[2]] == 0 ~ lv[2], TRUE ~ NA_character_))
  list(mat = mat, meta = meta_data, id_map = id_map, desc = desc_stats)
}

## Bias-estimation matrix (features x samples, neutral IDs A00001 ...)
bias_matrix <- function(bias_ps, samples) {
  m <- t(otu_matrix(prune_samples(samples, bias_ps)))[, samples, drop = FALSE]
  m <- m[rowSums(m) > 0, , drop = FALSE]
  rownames(m) <- sprintf("A%05d", seq_len(nrow(m)))
  m
}

## One ANCOM-BC2 fit. `bias` (optional): matrix of all features for the
## bias estimation; `pseudo`: pseudo-count added to all counts.
## Returns the estimable results (finite LFC and standard error) by feature ID.
fit_ancombc2 <- function(mat, meta_data, prv_cut, bias = NULL, pseudo = 0) {
  coef_name <- "groupTest"
  use_agg <- !is.null(bias) && "aggregate_data" %in% names(formals(ANCOMBC::ancombc2))
  O2 <- mat + pseudo
  O1 <- if (use_agg) bias + pseudo else NULL
  dropped <- character()
  notes <- character()
  for (iter in 1:20) {
    args <- list(
      data = if (use_agg) O1[setdiff(rownames(O1), dropped), , drop = FALSE]
             else O2[setdiff(rownames(O2), dropped), , drop = FALSE],
      taxa_are_rows = TRUE, meta_data = meta_data, tax_level = NULL,
      fix_formula = "group", rand_formula = NULL,
      p_adj_method = "holm", pseudo_sens = TRUE,
      prv_cut = prv_cut, lib_cut = 0, s0_perc = 0.05,
      group = NULL, struc_zero = FALSE, neg_lb = FALSE,
      alpha = cfg$ancom_alpha, n_cl = cfg$threads, verbose = FALSE,
      global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE)
    if (use_agg) args$aggregate_data <- O2
    set.seed(cfg$seed)
    fit <- tryCatch(withCallingHandlers(
      do.call(ANCOMBC::ancombc2, args),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (str_detect(msg, "NaN")) invokeRestart("muffleWarning")
        if (str_detect(msg, "number of taxa used for estimating")) {
          notes <<- c(notes, str_squish(msg)); invokeRestart("muffleWarning")
        }
      }),
      error = function(e) e)
    if (!inherits(fit, "error")) break
    msg <- conditionMessage(fit)
    bad <- if (str_detect(msg, "Zero variances")) {
      str_extract_all(str_extract(msg, "following taxa[\\s\\S]*$"), "[AF]\\d{5}")[[1]]
    } else character()
    bad <- setdiff(bad, dropped)
    if (length(bad) == 0) stop(fit)
    dropped <- c(dropped, bad)
    if (!use_agg && length(setdiff(rownames(O2), dropped)) < 2) stop(fit)
  }
  if (inherits(fit, "error")) stop(fit)
  res <- as_tibble(fit$res) |>
    transmute(id        = str_extract(taxon, "F\\d{5}"),
              lfc       = .data[[paste0("lfc_", coef_name)]],
              se        = .data[[paste0("se_", coef_name)]],
              W         = .data[[paste0("W_", coef_name)]],
              p         = .data[[paste0("p_", coef_name)]],
              passed_ss = .data[[paste0("passed_ss_", coef_name)]]) |>
    dplyr::filter(!is.na(id), is.finite(lfc), is.finite(se), se > 0)
  list(res = res, n_bias = if (use_agg) nrow(O1) - length(dropped) else nrow(O2),
       dropped = dropped, notes = unique(notes))
}

## ANCOM-BC2 on one phyloseq object and one condition;
## bias_ps: phyloseq object with all features for the bias estimation
run_ancombc2_safe <- function(physeq, label, ct, bias_ps = NULL) {
  cond <- ct$id
  prep <- prepare_da_input(physeq, ct)
  bias <- if (!is.null(bias_ps)) bias_matrix(bias_ps, colnames(prep$mat)) else NULL
  msgs <- character()

  ## 1. primary fit
  f1 <- tryCatch(fit_ancombc2(prep$mat, prep$meta, cfg$ancom_prv_cut, bias, pseudo = 0),
                 error = function(e) e)
  if (inherits(f1, "error")) {
    msgs <- c(msgs, paste0("primary fit: ", conditionMessage(f1)))
    log_msg("  ANCOM-BC2 (", cond, ", ", label, ") primary fit failed: ", conditionMessage(f1))
    f1 <- NULL
  }
  res1 <- if (is.null(f1)) NULL else mutate(f1$res, ancombc2_fit = "primary")

  ## 2. pseudo-count fit for every feature the primary fit could not estimate
  need <- setdiff(prep$desc$id, res1$id)
  res2 <- NULL; f2 <- NULL
  if (length(need) > 0) {
    f2 <- tryCatch(fit_ancombc2(prep$mat, prep$meta, prv_cut = 0, bias, pseudo = 1),
                   error = function(e) e)
    if (inherits(f2, "error")) {
      msgs <- c(msgs, paste0("pseudo-count fit: ", conditionMessage(f2)))
      log_msg("  ANCOM-BC2 (", cond, ", ", label, ") pseudo-count fit failed: ", conditionMessage(f2))
      f2 <- NULL
    } else {
      res2 <- f2$res |> dplyr::filter(id %in% need) |>
        mutate(ancombc2_fit = "pseudo-count 1")
      log_msg("  ANCOM-BC2 (", cond, ", ", label, "): ", nrow(res2), " of ", length(need),
              " features not estimable in the primary fit tested with pseudo-count 1")
    }
  }

  ## Combine; Holm adjustment over all tested features of the table
  res <- bind_rows(res1, res2)
  if (nrow(res) == 0) res <- tibble(id = character(), lfc = numeric(), se = numeric(),
                                    W = numeric(), p = numeric(), passed_ss = logical(),
                                    ancombc2_fit = character())
  res <- res |>
    mutate(q = p.adjust(p, method = "holm"),
           diff = q <= cfg$ancom_alpha,
           significant = diff & passed_ss)
  untested <- setdiff(prep$desc$id, res$id)
  if (length(untested) > 0) {
    log_msg("WARNING: ANCOM-BC2 (", cond, ", ", label, "): ", length(untested),
            " features could not be estimated by either fit (Wilcoxon test only)")
  }

  reason <- function(d) case_when(
    !is.na(d$absent_in) ~ paste0("zero counts in all ", d$absent_in, " samples"),
    !is.na(d$constant_in) ~ paste0("constant counts in ", d$constant_in),
    d$prevalence_all <= cfg$ancom_prv_cut ~
      paste0("present in ", 100 * cfg$ancom_prv_cut, "% of samples or fewer"),
    TRUE ~ "standard error not estimable without pseudo-count")
  all <- prep$desc |>
    left_join(res, by = "id") |>
    mutate(ancombc2_status = case_when(
             ancombc2_fit == "primary" ~ "Tested (primary fit)",
             ancombc2_fit == "pseudo-count 1" ~ paste0(
               "Tested with pseudo-count 1 (", reason(pick(everything())), ")"),
             TRUE ~ "Not estimable by ANCOM-BC2; Wilcoxon test only"),
           .after = taxon) |>
    select(-id) |>
    arrange(is.na(p), p, wilcox_p_nominal)

  fit <- list(res = dplyr::filter(all, !is.na(p)), all = all,
              n_input = nrow(prep$mat), n_samples = ncol(prep$mat))
  notes <- c(f1$notes, f2$notes)
  status <- tibble(
    Condition = ct$cond,
    Table = label,
    Comparison = ct$label,
    Features_input = nrow(prep$mat),
    Tested_primary = sum(all$ancombc2_fit %in% "primary"),
    Tested_pseudocount = sum(all$ancombc2_fit %in% "pseudo-count 1"),
    Not_tested = sum(is.na(all$p)),
    Bias_features = if (!is.null(f1)) f1$n_bias else if (!is.null(f2)) f2$n_bias else NA_integer_,
    Bias_features_removed_zero_variance = length(unique(c(f1$dropped, f2$dropped))),
    Samples = ncol(prep$mat),
    Status = if (length(untested) == 0) "OK - all features tested" else "Partly tested",
    Message = paste(c(msgs, notes), collapse = " | "))
  list(fit = fit, status = status)
}

## Summary of one condition's ANCOM-BC2 results per table (level / group)
ancom_summary_tbl <- function(all_tbl, key, ct) {
  tg <- lv_tag(ct$lv)
  all_tbl |>
    group_by(.data[[key]]) |>
    summarise(Features_input = n(),
              Features_tested = sum(!is.na(p)),
              Tested_primary = sum(ancombc2_fit %in% "primary"),
              Tested_pseudocount = sum(ancombc2_fit %in% "pseudo-count 1"),
              Not_tested = sum(is.na(p)),
              Wilcoxon_p05 = sum(wilcox_p_nominal < 0.05, na.rm = TRUE),
              Significant_q05 = sum(diff, na.rm = TRUE),
              Also_pass_sensitivity = sum(significant, na.rm = TRUE),
              !!paste0("Higher_in_", tg[2]) := sum(significant & lfc > 0, na.rm = TRUE),
              !!paste0("Higher_in_", tg[1]) := sum(significant & lfc < 0, na.rm = TRUE),
              .groups = "drop") |>
    mutate(Condition = ct$cond, Comparison = ct$label, .before = 1)
}

## Top-N features by ANCOM-BC2 p-value, regardless of significance; features
## not estimable by ANCOM-BC2 (if any) follow, ranked by Wilcoxon p
select_top <- function(all, n = cfg$top_n_da) {
  ancom <- all |> dplyr::filter(!is.na(p)) |> arrange(p, desc(abs(lfc)))
  wilc  <- all |> dplyr::filter(is.na(p), !is.na(wilcox_p_nominal)) |>
    arrange(wilcox_p_nominal)
  bind_rows(ancom, wilc) |>
    slice_head(n = n) |>
    mutate(rank = row_number(),
           test = if_else(!is.na(p), "ANCOM-BC2", "Wilcoxon"))
}

## Title suffix naming the test used to rank the selected features
top_by <- function(top) {
  if (all(top$test == "Wilcoxon")) " by Wilcoxon p-value"
  else if (any(top$test == "Wilcoxon")) " by ANCOM-BC2, then Wilcoxon p-value"
  else " by ANCOM-BC2 p-value"
}

## Boxplots of relative abundance (%) for selected features, one condition
plot_top_features <- function(rel_long, top, labels, title, ct) {
  cond <- ct$cond; lv <- ct$lv
  rel_long <- rel_long |>
    dplyr::filter(sample_id %in% ct$samples) |>
    mutate(!!cond := factor(as.character(.data[[cond]]), levels = lv))
  facet_lab <- setNames(
    if_else(top$test == "ANCOM-BC2",
            sprintf("%s\nLFC = %.2f%s\np = %s, q = %s", labels, top$lfc,
                    if_else(top$ancombc2_fit %in% "pseudo-count 1", " (pseudo 1)", ""),
                    signif(top$p, 2), signif(top$q, 2)),
            sprintf("%s\n%s\nWilcoxon p = %s%s", labels,
                    if_else(is.na(top$lfc), "ANCOM-BC2: n.e.",
                            sprintf("LFC = %.2f (SE n.e.)", top$lfc)),
                    signif(top$wilcox_p_nominal, 2),
                    if_else(is.na(top$absent_in), "",
                            paste0("\n(zero in all ", top$absent_in, ")")))),
    top$taxon)
  df <- rel_long |>
    dplyr::filter(taxon %in% top$taxon) |>
    mutate(panel = factor(facet_lab[taxon], levels = facet_lab))
  ggplot(df, aes(.data[[cond]], rel_abund, fill = .data[[cond]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
    geom_jitter(aes(colour = .data[[cond]]), width = 0.12, height = 0, size = cfg$point_size) +
    facet_wrap(~ panel, scales = "free_y", nrow = 1) +
    scale_x_discrete(labels = scales::label_wrap(12)) +
    scale_fill_manual(values = cond_cols[[cond]][lv]) +
    scale_colour_manual(values = cond_cols[[cond]][lv]) +
    labs(x = NULL, y = "Relative abundance (%)", title = title,
         caption = paste0("Top ", nrow(top), top_by(top), ", regardless of ",
                          "significance",
                          if (any(top$ancombc2_fit %in% "pseudo-count 1"))
                            paste0(".\npseudo 1: ANCOM-BC2 with a pseudo-count of 1, for features ",
                                   "not estimable without it (e.g. zero in all samples of ",
                                   "one group)") else "",
                          if (any(top$test == "Wilcoxon"))
                            paste0(".\nn.e.: not estimable by ANCOM-BC2 (zero variance within a ",
                                   "group, or standard error not estimable); for these, ",
                                   "two-sided Wilcoxon rank-sum p (nominal)") else "",
                          ".\nLFC: bias-corrected natural-log fold change, ", lv[2], " vs ", lv[1],
                          "; q: Holm-adjusted p")) +
    theme(legend.position = "none", strip.text = element_text(size = rel(0.8)))
}

## Per-sample relative abundance (%) of every taxon / feature
## (microViz::tax_transform + microViz::ps_melt)
level_abundance <- function(physeq) {
  physeq |>
    microViz::tax_transform("compositional", rank = "unique") |>
    microViz::ps_get() |>
    microViz::ps_melt() |>
    transmute(taxon = OTU, sample_id = Sample, rel_abund = 100 * Abundance) |>
    left_join(design, by = "sample_id")
}

## Relative abundance of the top-N taxa at one level in every sample, samples
## grouped by the condition (microViz::comp_barplot)
plot_relative_abundance_samples <- function(ps_level, level, cond, top_n = cfg$top_n_abund,
                                            title = NULL) {
  sample_order <- design$sample_id[design$sample_id %in% sample_names(ps_level)]
  microViz::comp_barplot(ps_level, tax_level = level, n_taxa = top_n,
                         sample_order = sample_order, facet_by = cond,
                         label = "SAMPLE", bar_width = 0.9, other_name = "Other",
                         merge_other = TRUE) +
    labs(x = "Sample", y = "Relative abundance", fill = level, title = title,
         caption = paste0("Top ", top_n, " taxa; relative abundance among taxa retained ",
                          "after filtering at this level")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

## Mean relative abundance per group of one condition at one level: relative
## abundances summed per group (phyloseq::merge_samples) and re-closed to 100%
## by comp_barplot, which equals the mean of the per-sample relative abundances
relative_abundance_panel <- function(ps_level, level, cond, top_n = cfg$top_n_abund) {
  ps_mean <- ps_level |>
    microViz::tax_transform("compositional", rank = "unique") |>
    microViz::ps_get() |>
    phyloseq::merge_samples(group = cond) |>
    suppressWarnings()
  microViz::comp_barplot(ps_mean, tax_level = level, n_taxa = top_n,
                         sample_order = intersect(conditions[[cond]], sample_names(ps_mean)),
                         label = "SAMPLE", bar_width = 0.7, other_name = "Other",
                         counts_warn = FALSE) +
    labs(x = NULL, y = "Mean relative abundance", fill = NULL, title = level) +
    theme(legend.text = element_text(size = rel(0.75)),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))   # group names never overlap
}

## ---- 4. Taxonomic tree (microViz taxatree) ---------------------------------

## Circular taxonomic tree of all taxa in a level object, from Domain at the
## centre to `level` at the tips, every node labelled with its taxon name.
## Node size is prevalence (proportion of samples in which the taxon occurs).
plot_taxonomy_tree <- function(ps_level, file, title = "") {
  strip_rank <- function(x) str_remove(x, "^[A-Za-z]: ")
  ps_t <- microViz::tax_prepend_ranks(ps_level)
  n_nodes <- sum(map_int(rank_names(ps_t), ~ n_distinct(as(tax_table(ps_t), "matrix")[, .x])))
  size_in <- max(10, min(30, 0.06 * n_nodes + 6))
  k <- size_in / 10
  p <- microViz::taxatree_plotkey(ps_t, rank != "Domain",
                                  .draw_label = FALSE, colour = "grey50",
                                  node_size_range = c(2, 6) * k, title = NULL) |>
    microViz::taxatree_plot_labels(taxon_renamer = strip_rank, size = 3.4 * k,
                                   fontface = "plain", seed = cfg$seed) +
    labs(title = title,
         caption = "Rings from centre: Domain > Phylum > ... ; node size = prevalence") +
    theme(plot.title = element_text(size = cfg$font_size * k, face = "bold"),
          plot.caption = element_text(size = (cfg$font_size - 3) * k))
  save_plot(p, dirname(file), tools::file_path_sans_ext(basename(file)), size_in, size_in,
            axes = FALSE, scale = FALSE)
  list(n_tips = ntaxa(ps_level), n_nodes = n_nodes)
}

## ---- 5. Taxonomic analysis -----------------------------------------------------

analyse_taxa <- function() {
  set.seed(cfg$seed)
  label <- "GTDB-Tk MAGs"
  root  <- step_dir(cfg$out_dir, "taxa")
  log_msg("==== Taxonomy: ", label, " ====")

  ## -- 01 Taxonomy: clean labels -----------------------------------------------
  d01 <- step_dir(root, "01_taxonomy")
  ps_raw <- round_counts(readRDS(ps_paths[["taxa"]]))
  log_msg("Imported ", ntaxa(ps_raw), " MAGs x ", nsamples(ps_raw), " samples; ",
          "taxonomy columns: ", paste(colnames(tax_table(ps_raw)), collapse = ", "),
          "; estimated reads rounded to whole reads")
  tx <- harmonise_taxonomy(ps_raw)
  save_table(tx$mapping,      d01, "taxonomy_column_mapping")
  save_table(tx$depth,        d01, "classification_depth")
  save_table(tx$label_log,    d01, "label_cleaning_by_rank")
  save_table(tx$placeholders, d01, "placeholder_labels_removed")
  ps_clean <- tx$ps

  ## -- 03 Filtering ------------------------------------------------------------
  ## (a) design samples, Bacteria only
  ## (b) duplicated taxon names removed rank by rank (validate_levels)
  ## (c) one table per level containing only valid, unique taxon names
  d03 <- step_dir(root, "03_filtering")
  ps_design <- attach_design(ps_clean, "ps_taxa.rds")
  filter_log <- tibble(Step = "After attaching metadata",
                       Features = ntaxa(ps_design), Reads = sum(sample_sums(ps_design)))
  ps_filt <- subset_taxa(ps_design, Domain == "Bacteria")
  ps_filt <- prune_taxa(taxa_sums(ps_filt) > 0, ps_filt)
  filter_log <- add_row(filter_log, Step = "Domain = Bacteria",
                        Features = ntaxa(ps_filt), Reads = sum(sample_sums(ps_filt))) |>
    mutate(Pct_reads_retained = round(100 * Reads / first(Reads), 2))
  save_table(filter_log, d03, "filtering_log_features")

  vl <- validate_levels(as.data.frame(as(tax_table(ps_filt), "matrix"),
                                      stringsAsFactors = FALSE))
  tax_valid <- vl$tax
  save_table(vl$duplicates, d03, "duplicated_taxa_removed")
  save_table(vl$parent_filled, d03, "unclassified_parents_filled_from_same_name")
  save_table(rownames_to_column(tax_valid, "feature"), d03, "taxonomy_table_validated")
  if (nrow(vl$duplicates) > 0) {
    log_msg("Duplicated taxon names removed: ",
            paste(names(table(vl$duplicates$Rank)),
                  map_int(split(vl$duplicates$name, vl$duplicates$Rank), n_distinct),
                  sep = " = ", collapse = ", "))
  }

  levels_all <- intersect(std_ranks[-1], names(tax_valid))
  lv <- map(set_names(levels_all), ~ build_level_ps(ps_filt, tax_valid, .x))
  level_log <- map_dfr(lv, "log")

  ## A level is analysed only if it has at least 2 valid taxa and at least 2
  ## samples in each group of every condition after filtering
  min_group_n <- map_int(lv, function(x) {
    sd <- data.frame(sample_data(x$ps))
    min(map_int(names(conditions), ~ min(table(factor(sd[[.x]], levels = conditions[[.x]])))))
  })
  level_log <- mutate(level_log, Min_samples_per_group = min_group_n[Level],
                      Analysed = Taxa >= 2 & Min_samples_per_group >= 2)
  save_table(level_log, d03, "filtering_log_by_level")
  too_small <- level_log$Level[!level_log$Analysed]
  if (length(too_small) > 0) {
    log_msg("Levels skipped (fewer than 2 valid taxa, or fewer than 2 samples in a ",
            "group after filtering): ", paste(too_small, collapse = ", "))
  }
  tax_levels <- setdiff(levels_all, too_small)
  ps_lv  <- map(lv[tax_levels], "ps")

  ## Taxonomy verification of every level object (used by every later step)
  tax_check <- imap(ps_lv, verify_level_taxonomy)
  ps_lv <- map(tax_check, "ps")
  tax_verif   <- map_dfr(tax_check, "report")
  tax_relabel <- map_dfr(tax_check, "relabel")
  save_table(tax_verif, d03, "taxonomy_name_verification")
  save_table(tax_relabel, d03, "taxonomy_relabelled")
  log_msg("Taxonomy names verified for ", length(ps_lv), " levels; ",
          nrow(tax_relabel), " labels re-labelled")

  ## Coverage audit: every taxon label in the design samples, analysed or not
  coverage <- taxa_coverage_audit(
    raw_labels   = tx$raw_labels,
    clean_tax    = as.data.frame(as(tax_table(ps_clean), "matrix"), stringsAsFactors = FALSE),
    design_feats = taxa_names(ps_design),
    bac_feats    = taxa_names(ps_filt),
    tax_valid    = tax_valid,
    level_taxa   = map(ps_lv, ~ as(tax_table(.x), "matrix")[, rank_names(.x)[ncol(tax_table(.x))]]))
  coverage_summary <- coverage |>
    group_by(Rank) |>
    summarise(Labels_in_design_samples = n(),
              Known_taxa = n_distinct(Name_analysed[Status == "Analysed"]) +
                sum(Status == "Known taxon - excluded"),
              Analysed = n_distinct(Name_analysed, na.rm = TRUE),
              Known_excluded_nonbacterial = sum(Status == "Known taxon - excluded" &
                                                  str_detect(Reason, "Not Bacteria")),
              Known_excluded_other = sum(Status == "Known taxon - excluded" &
                                           !str_detect(Reason, "Not Bacteria")),
              Not_known_labels = sum(Status == "Not a known taxon - excluded"),
              .groups = "drop") |>
    arrange(match(Rank, std_ranks))
  save_table(coverage, d03, "taxa_coverage_audit")
  save_table(coverage_summary, d03, "taxa_coverage_summary")
  log_msg("Coverage of known taxa (analysed / known): ",
          paste0(coverage_summary$Rank, " ", coverage_summary$Analysed, "/",
                 coverage_summary$Known_taxa, collapse = "; "))
  unexplained <- dplyr::filter(coverage, Status == "Known taxon - excluded",
                               !str_detect(Reason, "Not Bacteria"))
  if (nrow(unexplained) > 0) {
    log_msg("NOTE: known taxa not analysed for reasons other than non-bacterial: ",
            paste0(unexplained$Rank, ":", unexplained$Raw_label, " (", unexplained$Reason, ")",
                   collapse = "; "))
  }
  saveRDS(list(ps_filtered = ps_filt, tax_valid = tax_valid, levels = ps_lv),
          file.path(d03, "filtered_objects.rds"))

  ## Raw abundance tables (estimated reads, not rarefied or transformed), as TSV
  ## and as one Excel workbook:
  ##   features_unfiltered - every MAG of ps_taxa.rds, design samples, GTDB-Tk taxonomy
  ##   features_filtered   - bacterial MAGs, cleaned and validated taxonomy
  ##   <Level>             - estimated reads summed per valid taxon at each level
  try_step("Raw abundance export", {
    dra <- step_dir(d03, "raw_abundance")
    counts_tbl <- function(ps, tax) {
      cnt <- t(otu_matrix(ps))[, design$sample_id, drop = FALSE]     # features x samples
      as_tibble(tax[rownames(cnt), , drop = FALSE], rownames = "feature") |>
        bind_cols(as_tibble(cnt)) |>
        mutate(total_reads = rowSums(cnt))
    }
    raw_tbls <- list(
      features_unfiltered = counts_tbl(
        prune_samples(design$sample_id, ps_raw),
        as.data.frame(as(tax_table(ps_raw), "matrix"), stringsAsFactors = FALSE)),
      features_filtered = counts_tbl(ps_filt, tax_valid))
    lvl_tbls <- imap(lv, function(x, level) {
      cnt <- t(otu_matrix(x$ps))
      cnt <- cnt[, intersect(design$sample_id, colnames(cnt)), drop = FALSE]
      lin <- as.data.frame(as(tax_table(x$ps), "matrix"), stringsAsFactors = FALSE)
      as_tibble(lin[rownames(cnt), setdiff(names(lin), level), drop = FALSE]) |>
        mutate(!!level := rownames(cnt)) |>
        bind_cols(as_tibble(cnt)) |>
        mutate(total_reads = rowSums(cnt)) |>
        arrange(desc(total_reads))
    })
    raw_tbls <- c(raw_tbls, lvl_tbls)
    iwalk(raw_tbls, ~ write_tsv(.x, file.path(dra, paste0("raw_counts_", safe_name(.y), ".tsv")),
                                na = ""))
    save_excel(raw_tbls, dra, "raw_counts_gtdbtk")
    log_msg("Raw abundance tables written: ", paste(names(raw_tbls), collapse = ", "))
  })
  log_msg("Per-level tables: ",
          paste0(level_log$Level, " ", level_log$Taxa, " taxa (",
                 level_log$Pct_reads_retained, "% reads)", collapse = "; "))

  ## -- 04 Library sizes (overall: all bacterial MAGs) --------------------------
  d04 <- step_dir(root, "04_library_size")
  lib <- tibble(sample_id = sample_names(ps_filt), reads = sample_sums(ps_filt)) |>
    left_join(design, by = "sample_id") |>
    mutate(sample_lab = factor(sample_id, levels = design$sample_id))
  lib_groups <- map_dfr(names(conditions), function(cond) {
    lib |>
      group_by(Group = .data[[cond]]) |>
      summarise(n = n(), Min = min(reads), Median = median(reads), Max = max(reads),
                .groups = "drop") |>
      mutate(Condition = cond, Group = as.character(Group), .before = 1)
  })
  ## two-sided Wilcoxon rank-sum test for every pair of levels (nominal p)
  lib_tests <- map_dfr(contrasts, function(ct) {
    x <- lib$reads[lib[[ct$cond]] == ct$lv[1]]; y <- lib$reads[lib[[ct$cond]] == ct$lv[2]]
    tibble(Condition = ct$cond, group1 = ct$lv[1], group2 = ct$lv[2],
           n1 = length(x), n2 = length(y),
           Wilcoxon_p = wilcox.test(x, y, exact = FALSE)$p.value)
  })
  save_table(select(lib, -sample_lab), d04, "library_sizes")
  save_table(lib_groups, d04, "library_sizes_by_condition")
  save_table(lib_tests, d04, "library_sizes_wilcoxon")
  for (cond in names(conditions)) {
    try_step(paste0("Library size plot (", cond, ")"), {
      lt <- dplyr::filter(lib_tests, Condition == cond)
      ## samples grouped by the condition (reference level first), then in
      ## sample order within each group
      lib_c <- mutate(lib, sample_lab = factor(sample_id, levels = sample_id[
        order(as.integer(.data[[cond]]), match(sample_id, design$sample_id))]))
      p_lib <- ggplot(lib_c, aes(sample_lab, reads, fill = .data[[cond]])) +
        geom_col() +
        scale_fill_manual(values = cond_cols[[cond]]) +
        scale_y_continuous(labels = scales::comma) +
        labs(x = "Sample", y = "Estimated reads (bacterial MAGs)", fill = cond_title(cond),
             title = label,
             subtitle = paste0("Wilcoxon (nominal p): ",
                               paste(sprintf("%s vs %s: p = %.3g", lt$group1, lt$group2,
                                             lt$Wilcoxon_p), collapse = "\n"))) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      save_plot(p_lib, d04, paste0("library_sizes_", cond), 10, 5.5)
    })
  }

  ## -- 05 Alpha rarefaction (overall) ------------------------------------------
  d05 <- step_dir(root, "05_alpha_rarefaction")
  otu_mat <- otu_matrix(ps_filt)                   # samples x MAGs, plain matrix
  rare_depth <- min(sample_sums(ps_filt))
  try_step("Rarefaction curves", {
    rc <- rarefaction_curves(otu_mat)
    save_table(rc, d05, "rarefaction_curves")
    save_table(attr(rc, "saturation"), d05, "rarefaction_saturation_depth")
    for (cond in names(conditions)) {
      save_plot(plot_rarefaction(rc, cond, rare_depth, label, "Estimated reads sampled"),
                d05, paste0("alpha_rarefaction_curves_", cond), 9.5, 6.5)
    }
  })
  ps_rare <- rarefy_ps(ps_filt, rare_depth)
  saveRDS(ps_rare, file.path(d05, "ps_rarefied.rds"))
  save_table(tibble(Object = c("Before rarefying", "After rarefying"),
                    Features = c(ntaxa(ps_filt), ntaxa(ps_rare)),
                    Reads_per_sample = c(paste(range(sample_sums(ps_filt)), collapse = "-"),
                                         as.character(rare_depth))),
             d05, "rarefaction_effect")
  log_msg("Rarefied to ", format(rare_depth, big.mark = ","), " reads; ",
          ntaxa(ps_rare), " MAGs retained")

  ## -- 06 Alpha diversity (overall; Wilcoxon, nominal p) -----------------------
  d06 <- step_dir(root, "06_alpha_diversity")
  a_val <- alpha_values(ps_rare)
  save_table(a_val, d06, "alpha_diversity")
  al_by <- map(set_names(names(conditions)), ~ alpha_tests(a_val, .x))
  alpha_tests_all <- map_dfr(al_by, ~ mutate(.x$tests, index = as.character(index)))
  save_table(map_dfr(al_by, "summary") |> mutate(index = as.character(index)),
             d06, "alpha_diversity_summary")
  save_table(alpha_tests_all, d06, "alpha_wilcoxon")
  iwalk(al_by, function(al, cond) {
    try_step(paste0("Alpha diversity plot (", cond, ")"),
             save_plot(plot_alpha(al, cond, paste0(label, " - ", cond_title(cond))),
                       d06, paste0("alpha_diversity_", cond), 12, 8))
  })
  log_msg("Alpha diversity: ", sum(alpha_tests_all$p < 0.05, na.rm = TRUE), " of ",
          sum(!is.na(alpha_tests_all$p)), " index x condition tests with nominal p < 0.05")

  ## -- 07 Relative abundance (tables per level; figures per condition) --------
  d07 <- step_dir(root, "07_relative_abundance")
  rel_ab <- map(ps_lv, level_abundance)
  iwalk(rel_ab, function(ab, level) {
    save_table(ab |>
                 select(taxon, sample_id, rel_abund) |>
                 pivot_wider(names_from = sample_id, values_from = rel_abund) |>
                 mutate(mean_all = rowMeans(across(-taxon))) |>
                 arrange(desc(mean_all)),
               d07, paste0("relative_abundance_", safe_name(level)))
  })
  for (cond in names(conditions)) {
    d07c <- step_dir(d07, cond)
    iwalk(ps_lv, function(ps_l, level) {
      try_step(paste0("Relative abundance per sample (", cond, ", ", level, ")"),
               save_plot(plot_relative_abundance_samples(
                           ps_l, level, cond, title = paste0(level, " - ", label)),
                         d07c, paste0("relative_abundance_samples_", safe_name(level)), 15, 9))
    })
    try_step(paste0("Relative abundance overview plot (", cond, ")"), {
      panels <- imap(ps_lv, ~ relative_abundance_panel(.x, .y, cond))
      p_rel <- wrap_plots(panels, ncol = 3) +
        plot_annotation(
          title = paste0("Mean relative abundance by ", tolower(cond_title(cond)), " - ", label),
          caption = paste0("Top ", cfg$top_n_abund, " taxa per level; relative abundance ",
                           "among taxa retained after filtering at each level"))
      save_plot(p_rel, d07c, "relative_abundance_overview", 21,
                6.5 * ceiling(length(panels) / 3) + 1, k = 1.5)
    })
  }

  ## -- 08 Beta diversity (overall; Bray-Curtis PCoA) ---------------------------
  d08 <- step_dir(root, "08_beta_diversity")
  try_step("PCoA", {
    b <- pcoa_bray(ps_filt)
    save_table(b$coords, d08, "pcoa_bray_coordinates")
    saveRDS(b$dist, file.path(d08, "bray_curtis_distance.rds"))
    for (cond in names(conditions)) {
      p <- plot_pcoa(b, cond, paste0("Bray-Curtis PCoA - ", label))
      if (!is.null(p)) save_plot(p, d08, paste0("pcoa_bray_", cond), 8.5, 7)
    }
  })

  ## -- 09 PERMANOVA + PERMDISP (all features and each level, per condition) ----
  d09 <- step_dir(root, "09_permanova")
  perm_levels <- c("All features", tax_levels)
  ## all levels of each condition together, plus every pair of levels for
  ## conditions with more than two levels
  perm_runs <- c(map(names(conditions), ~ list(cond = .x, ct = NULL)),
                 map(keep(contrasts, ~ length(conditions[[.x$cond]]) > 2),
                     ~ list(cond = .x$cond, ct = .x)))
  permanova <- map_dfr(perm_runs, function(r) {
    d09r <- step_dir(d09, "raw_output", r$cond)
    imap_dfr(c(list(`All features` = ps_filt), ps_lv),
             ~ try_step(paste0("PERMANOVA (", r$cond, if (!is.null(r$ct)) paste0(", ", r$ct$label),
                               ", ", .y, ")"),
                        permanova_bray(.x, .y, r$cond, raw_dir = d09r, ct = r$ct)))
  })
  if (nrow(permanova) == 0) {
    permanova <- tibble(Condition = character(), Comparison = character(), Level = character(),
                        R2 = numeric(), p = numeric(), PERMDISP_p = numeric())
  }
  permanova <- permanova |>
    group_by(Condition, Comparison) |>
    mutate(p_adj = p.adjust(p, method = "BH"),
           PERMDISP_p_adj = p.adjust(PERMDISP_p, method = "BH")) |>
    ungroup() |>
    mutate(Level = factor(Level, levels = perm_levels))
  save_table(permanova, d09, "permanova_by_level")
  try_step("PERMANOVA plot", {
    p_perm <- ggplot(dplyr::filter(permanova, Level != "All features",
                                   Comparison == "All groups") |> droplevels(),
                     aes(Level, R2)) +
      geom_col(fill = "grey60") +
      geom_text(aes(label = paste0("p = ", signif(p, 2))), vjust = -0.4, size = cfg$text_size) +
      facet_wrap(~ Condition, ncol = 1, labeller = as_labeller(cond_title)) +
      scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
      labs(x = NULL, y = expression(R^2 ~ "(variation explained by the condition)"), title = label,
           caption = paste0("PERMANOVA on Bray-Curtis, ", cfg$n_perm, " permutations")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    save_plot(p_perm, d09, "permanova_R2_by_level", 9, 11)
  })

  ## -- 10 ANCOM-BC2 (per level, per condition) ---------------------------------
  d10 <- step_dir(root, "10_ancombc2")
  save_table(tax_verif, d10, "taxonomy_name_verification")
  save_table(tax_relabel, d10, "taxonomy_relabelled")
  ancom_ct <- map(contrasts, function(ct) {
    cond <- ct$cond
    d10c <- step_dir(d10, ct$id)
    runs <- map(set_names(tax_levels), function(level) {
      log_msg("ANCOM-BC2 (", cond, ", ", ct$label, "): ", level, " (", ntaxa(ps_lv[[level]]),
              " taxa, ", length(ct$samples), " samples)")
      run_ancombc2_safe(ps_lv[[level]], level, ct, bias_ps = ps_filt)
    })
    ancom <- compact(map(runs, "fit"))
    status <- bind_rows(
      map_dfr(runs, "status"),
      level_log |>
        dplyr::filter(!Analysed) |>
        transmute(Condition = cond, Comparison = ct$label, Table = Level, Features_input = Taxa,
                  Samples = Samples, Status = "Skipped",
                  Message = "Fewer than 2 valid taxa or fewer than 2 samples per group")) |>
      rename(Level = Table) |>
      arrange(match(Level, std_ranks))
    save_table(status, d10c, "ancombc2_run_status")
    all_tbl <- imap_dfr(ancom, ~ mutate(.x$all, Level = .y, .before = 1))
    save_table(all_tbl, d10c, "ancombc2_by_level")
    summ <- if (nrow(all_tbl) > 0) {
      ancom_summary_tbl(all_tbl, "Level", ct) |> arrange(match(Level, std_ranks))
    } else tibble()
    save_table(summ, d10c, "ancombc2_summary")

    ## Excel workbook: Summary + Run status + one sheet per level, each with the
    ## higher-rank lineage, prevalence and abundance per group and all statistics
    try_step(paste0("ANCOM-BC2 Excel export (", ct$id, ")"), {
      level_sheets <- imap(ancom, function(a, level) {
        lineage <- as.data.frame(as(tax_table(ps_lv[[level]]), "matrix"),
                                 stringsAsFactors = FALSE) |>
          rownames_to_column("taxon")
        lineage <- lineage[, setdiff(names(lineage), level), drop = FALSE]
        a$all |>
          left_join(lineage, by = "taxon") |>
          relocate(any_of(std_ranks), .after = taxon) |>
          rename(!!level := taxon)
      })
      save_excel(c(list(Summary = summ, Run_status = status), level_sheets),
                 d10c, "ancombc2_results_by_level")
    })

    ## Top-N taxa per level by p-value, regardless of significance
    d10t <- step_dir(d10c, paste0("top", cfg$top_n_da, "_taxa"))
    top_da <- imap_dfr(ancom, ~ mutate(select_top(.x$all), Level = .y, .before = 1))
    save_table(top_da, d10t, paste0("top", cfg$top_n_da, "_taxa_by_level"))
    for (level in unique(top_da$Level)) {
      try_step(paste0("Top-", cfg$top_n_da, " plot (", ct$id, ", ", level, ")"), {
        top <- dplyr::filter(top_da, Level == level)
        p_top <- plot_top_features(rel_ab[[level]], top, str_wrap(top$taxon, 18),
                                   paste0("Top ", nrow(top), " ", level, top_by(top), " - ",
                                          cond_title(cond), ": ", ct$label), ct)
        save_plot(p_top, d10t, paste0("top", cfg$top_n_da, "_", safe_name(level)),
                  max(6, 4.2 * nrow(top) + 1), 7)
      })
    }
    list(all = all_tbl, summary = summ, top = top_da)
  })

  ## -- 11 Taxonomic tree -------------------------------------------------------
  d11 <- step_dir(root, "11_taxonomic_tree")
  ## One tree per lowest rank: Domain -> Genus (all valid genera) and
  ## Domain -> Species (all valid species)
  tree_ranks <- intersect(c("Genus", "Species"), tax_levels)
  tree_res <- NULL
  for (tr_rank in tree_ranks) {
    tr <- try_step(paste0("Taxonomic tree (", tr_rank, ")"), {
      res_tr <- plot_taxonomy_tree(
        ps_lv[[tr_rank]],
        file.path(d11, paste0("taxonomic_tree_domain_to_", safe_name(tr_rank), ".png")),
        title = paste0("All identified taxa, Domain to ", tr_rank, " - ", label))
      save_table(as.data.frame(as(tax_table(ps_lv[[tr_rank]]), "matrix")) |>
                   arrange(across(everything())),
                 d11, paste0("taxonomic_lineages_to_", safe_name(tr_rank)))
      log_msg("Taxonomic tree to ", tr_rank, ": ", res_tr$n_tips, " tips, ",
              res_tr$n_nodes, " nodes")
      res_tr
    })
    if (tr_rank == "Genus" || is.null(tree_res)) tree_res <- tr
  }

  ## Heat trees (metacoder): composition, and per condition the differential
  ## tree (Wilcoxon, nominal p)
  if (has_pkg[["metacoder"]] && "Genus" %in% names(ps_lv)) {
    tryCatch({
      ## metacoder's rank functions look up its dataset 'ranks_ref' by name,
      ## visible only when metacoder is attached; loaded explicitly here
      if (!exists("ranks_ref", envir = globalenv(), inherits = FALSE)) {
        utils::data("ranks_ref", package = "metacoder", envir = globalenv())
      }
      obj <- metacoder::parse_phyloseq(ps_lv[["Genus"]])
      sample_ids <- obj$data$sample_data$sample_id
      obj$data$tax_prop  <- metacoder::calc_obs_props(obj, "otu_table", cols = sample_ids)
      obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "tax_prop", cols = sample_ids)
      obj$data$tax_summary <- tibble(
        taxon_id  = obj$data$tax_abund$taxon_id,
        mean_prop = rowMeans(obj$data$tax_abund[, sample_ids]))

      set.seed(cfg$seed)
      p_ht1 <- metacoder::heat_tree(
        obj, node_label = taxon_names, node_size = n_obs,
        node_color = mean_prop * 100,
        node_color_range = c("grey90", "#FDB863", "#E66101", "#5E3C99"),
        node_size_axis_label = "Number of genera",
        node_color_axis_label = "Mean relative abundance (%)",
        node_label_size_range = c(0.012, 0.025),
        layout = "davidson-harel", initial_layout = "reingold-tilford")
      save_plot(p_ht1, d11, "heat_tree_composition", 14, 14, axes = FALSE, scale = FALSE)

      for (ct in contrasts) {
        cond <- ct$cond
        ct_cols <- sample_ids[sample_ids %in% ct$samples]
        groups <- as.character(obj$data$sample_data[[cond]][match(ct_cols, sample_ids)])
        obj$data$diff_table <- metacoder::compare_groups(obj, data = "tax_abund",
                                                         cols = ct_cols, groups = groups)
        obj$data$diff_table <- obj$data$diff_table |>
          mutate(log2_ratio_plot = case_when(
                   is.na(wilcox_p_value) | wilcox_p_value >= 0.05 ~ 0,
                   is.infinite(log2_median_ratio) ~ sign(log2_median_ratio) * 5,
                   TRUE ~ pmax(pmin(log2_median_ratio, 5), -5)))
        contrast <- distinct(obj$data$diff_table, treatment_1, treatment_2)
        save_table(left_join(obj$data$diff_table,
                             tibble(taxon_id = names(obj$taxon_names()),
                                    taxon = obj$taxon_names()), by = "taxon_id"),
                   d11, paste0("heat_tree_wilcoxon_", ct$id))
        set.seed(cfg$seed)
        p_ht2 <- metacoder::heat_tree(
          obj, node_label = taxon_names, node_size = n_obs,
          node_color = log2_ratio_plot, node_color_interval = c(-5, 5),
          node_color_range = c(cond_cols[[cond]][[contrast$treatment_2]], "grey85",
                               cond_cols[[cond]][[contrast$treatment_1]]),
          node_size_axis_label = "Number of genera",
          node_color_axis_label = paste0("log2 median ratio (", contrast$treatment_1,
                                         " / ", contrast$treatment_2, "), Wilcoxon p < 0.05"),
          node_label_size_range = c(0.012, 0.025),
          layout = "davidson-harel", initial_layout = "reingold-tilford")
        save_plot(p_ht2, d11, paste0("heat_tree_differential_", ct$id), 14, 14,
                  axes = FALSE, scale = FALSE)
      }
      log_msg("Heat trees written")
    }, error = function(e) log_msg("Heat trees failed: ", conditionMessage(e)))
  } else if (!has_pkg[["metacoder"]]) {
    log_msg("metacoder not installed; heat trees skipped")
  }

  ## -- Summary -----------------------------------------------------------------
  ## PERMANOVA row of one comparison: all groups for a two-level condition,
  ## the pairwise test otherwise
  perm_cmp <- function(ct) if (length(conditions[[ct$cond]]) == 2) "All groups" else ct$label
  summary_overall <- map_dfr(contrasts, function(ct) {
    cond <- ct$cond
    is_pair <- function(t) t$Condition == cond & t$group1 == ct$lv[1] & t$group2 == ct$lv[2]
    pa <- dplyr::filter(permanova, Condition == cond, Level == "All features")
    tibble(
      Condition = cond,
      Comparison = ct$label,
      Samples = length(ct$samples),
      MAGs = ntaxa(ps_filt),
      Rarefaction_depth = rare_depth,
      Library_size_wilcox_p = lib_tests$Wilcoxon_p[is_pair(lib_tests)][1],
      Alpha_indices_p05 = sum(alpha_tests_all$p[is_pair(alpha_tests_all)] < 0.05, na.rm = TRUE),
      PERMANOVA_R2_all_features = pa$R2[pa$Comparison == perm_cmp(ct)][1],
      PERMANOVA_p_all_features  = pa$p[pa$Comparison == perm_cmp(ct)][1],
      PERMANOVA_all_groups_R2   = pa$R2[pa$Comparison == "All groups"][1],
      PERMANOVA_all_groups_p    = pa$p[pa$Comparison == "All groups"][1],
      Tree_tips = if (is.null(tree_res)) NA_integer_ else tree_res$n_tips)
  })
  save_table(summary_overall, root, "summary_overall")

  summary_tbl <- try_step("Summary table", map_dfr(contrasts, function(ct) {
    cond <- ct$cond
    level_log |>
      dplyr::filter(Level %in% tax_levels) |>
      select(Level, Taxa, Pct_reads_retained, Samples) |>
      left_join(permanova |> dplyr::filter(Condition == cond, Comparison == perm_cmp(ct)) |>
                  mutate(Level = as.character(Level)) |>
                  select(Level, PERMANOVA_R2 = R2, PERMANOVA_p = p), by = "Level") |>
      left_join(ancom_ct[[ct$id]]$summary |>
                  select(any_of(c("Level", "Features_tested", "Also_pass_sensitivity"))) |>
                  rename(any_of(c(ANCOMBC2_tested = "Features_tested",
                                  ANCOMBC2_significant = "Also_pass_sensitivity"))),
                by = "Level") |>
      mutate(Condition = cond, Comparison = ct$label, .before = 1)
  }))
  if (is.null(summary_tbl)) summary_tbl <- tibble()
  save_table(summary_tbl, root, "summary_by_level")

  results <- list(tax_valid = tax_valid, ps_filtered = ps_filt, ps_rarefied = ps_rare,
                  levels = ps_lv, classification_depth = tx$depth,
                  summary_overall = summary_overall, alpha_tests = alpha_tests_all,
                  permanova = permanova, ancombc2 = ancom_ct, summary = summary_tbl)
  saveRDS(results, file.path(root, "results.rds"))
  results
}

## ---- 6. Functional analysis --------------------------------------------------

## Filtering of one functional table (as for PICRUSt2 tables in the 16S
## workflow). Excluded, in this order:
##   1. rows with a missing feature ID or name
##   2. names matching invalid_function_regex (unknown/unclassified/unnamed)
##   3. duplicated feature IDs (all copies)
##   4. names shared by more than one feature ID (all copies;
##      cfg$drop_duplicated_descriptions)
##   5. features with zero abundance in every design sample
## Abundances are estimated reads, rounded to whole reads.
filter_functions <- function(ps_g, group) {
  tt  <- as.data.frame(as(tax_table(ps_g), "matrix"), stringsAsFactors = FALSE)
  if (!all(c("ID", "NAME") %in% names(tt))) stop(group, ": tax_table needs columns ID and NAME")
  cnt <- t(otu_matrix(ps_g))                                  # features x samples
  samples <- design$sample_id
  miss <- setdiff(samples, colnames(cnt))
  if (length(miss) > 0) stop(group, " table lacks samples: ", paste(miss, collapse = ", "))

  df <- tibble(feature = str_trim(tt[rownames(cnt), "ID"]),
               description = str_trim(tt[rownames(cnt), "NAME"])) |>
    bind_cols(as_tibble(cnt[, samples, drop = FALSE])) |>
    mutate(total = rowSums(across(all_of(samples))))
  total_all <- sum(df$total)
  invalid_re <- regex(cfg$invalid_function_regex, ignore_case = TRUE)

  df <- df |>
    mutate(reason = case_when(
      is.na(feature) | feature == "" | is.na(description) | description == "" ~
        "1 Missing ID or name",
      str_detect(description, invalid_re) ~ "2 Unknown/unclassified/unnamed name",
      feature %in% feature[duplicated(feature)] ~ "3 Duplicated feature ID",
      TRUE ~ NA_character_))
  valid_desc <- df$description[is.na(df$reason)]
  dup_desc   <- if (cfg$drop_duplicated_descriptions) unique(valid_desc[duplicated(valid_desc)]) else character()
  df <- df |>
    mutate(reason = case_when(
      !is.na(reason) ~ reason,
      description %in% dup_desc ~ "4 Duplicated name",
      total == 0 ~ "5 Zero abundance in all samples",
      TRUE ~ "Retained"))

  log <- df |>
    group_by(reason) |>
    summarise(Features = n(), Abundance = sum(total), .groups = "drop") |>
    mutate(Group = group, Pct_abundance = round(100 * Abundance / total_all, 2),
           .before = 1)
  kept <- dplyr::filter(df, reason == "Retained")
  mat  <- as.matrix(kept[, samples])
  rownames(mat) <- kept$feature
  mat <- round(mat)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  list(mat = mat,
       desc = setNames(kept$description, kept$feature),
       log = log,
       excluded = dplyr::filter(df, reason != "Retained") |>
         select(feature, description, reason, total))
}

analyse_functions <- function() {
  ps_function <- readRDS(ps_paths[["function_groups"]])
  miss_g <- setdiff(fun_groups, names(ps_function))
  if (length(miss_g) > 0) {
    stop("Functional groups missing from ps_function.rds: ", paste(miss_g, collapse = ", "),
         "; present: ", paste(names(ps_function), collapse = ", "))
  }
  root <- step_dir(cfg$out_dir, "function")
  log_msg("==== Functional analysis: ", paste(fun_groups, collapse = ", "), " ====")

  ## -- 02 Import and filtering -------------------------------------------------
  d01 <- step_dir(root, "02_filtering")
  fl <- map(set_names(fun_groups), function(group) {
    ps_g <- attach_design(ps_function[[group]], paste0("ps_function.rds$", group))
    f <- filter_functions(ps_g, group)
    save_table(f$excluded, d01, paste0("excluded_features_", safe_name(group)))
    save_table(as.data.frame(f$mat) |>
                 rownames_to_column("feature") |>
                 mutate(description = f$desc[feature], .after = feature),
               d01, paste0("filtered_abundance_", safe_name(group)))
    log_msg(group, ": ", nrow(f$mat), " features retained of ", sum(f$log$Features))
    f
  })
  save_table(map_dfr(fl, "log"), d01, "filtering_log")

  ## phyloseq containers (features as "taxa") for the shared helpers
  ps_fun <- map(fl, function(f) {
    tt <- matrix(rownames(f$mat), ncol = 1, dimnames = list(rownames(f$mat), "Function"))
    phyloseq(otu_table(f$mat, taxa_are_rows = TRUE), tax_table(tt), sam_data)
  })

  ## -- 03 Alpha diversity (Wilcoxon, nominal p) --------------------------------
  ## Computed on the filtered functional abundances, not rarefied (as for the
  ## PICRUSt2 tables in the 16S workflow)
  d03 <- step_dir(root, "03_alpha_diversity")
  fun_alpha <- imap(ps_fun, function(ps_g, group) {
    try_step(paste0("Functional alpha diversity (", group, ")"), {
      a_val <- alpha_values(ps_g)
      save_table(a_val, d03, paste0("alpha_diversity_", safe_name(group)))
      al_by <- map(set_names(names(conditions)), ~ alpha_tests(a_val, .x))
      save_table(map_dfr(al_by, "summary") |> mutate(index = as.character(index)),
                 d03, paste0("alpha_summary_", safe_name(group)))
      iwalk(al_by, function(al, cond) {
        save_plot(plot_alpha(al, cond, paste0("Functional alpha diversity - ", group, " - ",
                                              cond_title(cond))),
                  d03, paste0("alpha_diversity_", safe_name(group), "_", cond), 12, 8)
      })
      map_dfr(al_by, ~ mutate(.x$tests, index = as.character(index)))
    })
  }) |> compact()
  save_table(imap_dfr(fun_alpha, ~ mutate(.x, Group = .y, .before = 1)),
             d03, "alpha_wilcoxon_all_groups")

  ## -- 04 Beta diversity (Bray-Curtis PCoA; PCA as ggpicrust2::pathway_pca) ----
  d04 <- step_dir(root, "04_beta_diversity")
  iwalk(ps_fun, function(ps_g, group) {
    try_step(paste0("Functional PCoA (", group, ")"), {
      b <- pcoa_bray(ps_g)
      save_table(b$coords, d04, paste0("pcoa_bray_", safe_name(group)))
      for (cond in names(conditions)) {
        p <- plot_pcoa(b, cond, paste0("Bray-Curtis PCoA - ", group))
        if (!is.null(p)) save_plot(p, d04, paste0("pcoa_bray_", safe_name(group), "_", cond), 8.5, 7)
      }
    })
    for (cond in names(conditions)) {
      try_step(paste0("Functional PCA (", group, ", ", cond, ")"),
               plot_pca(ps_g, cond, d04, paste0("pca_", safe_name(group), "_", cond)))
    }
  })

  ## -- 05 PERMANOVA + PERMDISP -------------------------------------------------
  d05 <- step_dir(root, "05_permanova")
  ## all levels of each condition together, plus every pair of levels for
  ## conditions with more than two levels
  perm_runs <- c(map(names(conditions), ~ list(cond = .x, ct = NULL)),
                 map(keep(contrasts, ~ length(conditions[[.x$cond]]) > 2),
                     ~ list(cond = .x$cond, ct = .x)))
  fun_perm <- map_dfr(perm_runs, function(r) {
    d05r <- step_dir(d05, "raw_output", r$cond)
    imap_dfr(ps_fun, ~ try_step(paste0("Functional PERMANOVA (", r$cond,
                                       if (!is.null(r$ct)) paste0(", ", r$ct$label), ", ", .y, ")"),
                                permanova_bray(.x, .y, r$cond, raw_dir = d05r, ct = r$ct)))
  }) |>
    rename(Group = Level) |>
    group_by(Condition, Comparison) |>
    mutate(p_adj = p.adjust(p, method = "BH"),
           PERMDISP_p_adj = p.adjust(PERMDISP_p, method = "BH")) |>
    ungroup() |>
    mutate(Group = factor(Group, levels = fun_groups))
  save_table(fun_perm, d05, "permanova_by_group")
  try_step("Functional PERMANOVA plot",
    save_plot(ggplot(dplyr::filter(fun_perm, Comparison == "All groups"), aes(Group, R2)) +
                geom_col(fill = "grey60") +
                geom_text(aes(label = paste0("p = ", signif(p, 2))), vjust = -0.4,
                          size = cfg$text_size) +
                facet_wrap(~ Condition, ncol = 1, labeller = as_labeller(cond_title)) +
                scale_y_continuous(labels = scales::percent,
                                   expand = expansion(mult = c(0, 0.15))) +
                labs(x = NULL, y = expression(R^2 ~ "(variation explained by the condition)"),
                     title = "PERMANOVA on functional profiles",
                     caption = paste0("Bray-Curtis, ", cfg$n_perm, " permutations")) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
              d05, "permanova_R2_by_group", 9, 11))

  ## -- 06 ANCOM-BC2 and top-N features (per condition) -------------------------
  d06 <- step_dir(root, "06_ancombc2")
  fun_rel <- map(ps_fun, level_abundance)
  fun_ct <- map(contrasts, function(ct) {
    cond <- ct$cond
    d06c <- step_dir(d06, ct$id)
    d06t <- step_dir(d06c, paste0("top", cfg$top_n_da, "_features"))
    runs <- imap(fl, function(f, group) {
      log_msg("ANCOM-BC2 (", cond, ", ", ct$label, "): ", group, " (", nrow(f$mat), " features)")
      ## KEGG Module / Pathway abundances come from the KO annotations, so the
      ## KO features are used for their bias estimation; the other groups have
      ## enough features of their own
      bias <- if (group %in% c("KEGG_Module", "KEGG_Pathway")) ps_fun$KEGG_KO else NULL
      run_ancombc2_safe(ps_fun[[group]], group, ct, bias_ps = bias)
    })
    ancom  <- compact(map(runs, "fit"))
    status <- map_dfr(runs, "status") |> rename(Group = Table)
    save_table(status, d06c, "ancombc2_run_status")

    res <- imap_dfr(ancom, function(a, group) {
      a$all |>
        rename(feature = taxon) |>
        mutate(Group = group, description = fl[[group]]$desc[feature], .before = 1)
    })
    save_table(res, d06c, "ancombc2_by_group")
    summ <- if (nrow(res) > 0) ancom_summary_tbl(res, "Group", ct) else tibble()
    save_table(summ, d06c, "ancombc2_summary")

    try_step(paste0("Functional ANCOM-BC2 Excel export (", ct$id, ")"), {
      group_sheets <- imap(ancom, function(a, group) {
        a$all |>
          rename(feature = taxon) |>
          mutate(description = fl[[group]]$desc[feature], .after = feature)
      })
      save_excel(c(list(Summary = summ, Run_status = status), group_sheets), d06c,
                 "ancombc2_results_by_group")
    })

    top_fun <- imap_dfr(ancom, function(a, group) {
      select_top(a$all) |>
        mutate(Group = group, description = fl[[group]]$desc[taxon], .before = 1)
    })
    save_table(rename(top_fun, feature = taxon), d06t,
               paste0("top", cfg$top_n_da, "_features_by_group"))
    for (group in unique(top_fun$Group)) {
      try_step(paste0("Functional top-", cfg$top_n_da, " plot (", ct$id, ", ", group, ")"), {
        top <- dplyr::filter(top_fun, Group == group)
        labels <- paste0(top$taxon, "\n", str_wrap(str_trunc(top$description, 70), 24))
        p_top <- plot_top_features(fun_rel[[group]], top, labels,
                                   paste0("Top ", nrow(top), " ", group, " features", top_by(top),
                                          " - ", cond_title(cond), ": ", ct$label), ct)
        save_plot(p_top, d06t, paste0("top", cfg$top_n_da, "_", safe_name(group)),
                  max(6, 4.4 * nrow(top) + 1), 8)
      })
    }
    list(all = res, summary = summ, top = top_fun, status = status)
  })

  saveRDS(list(filtered = map(fl, ~ .x[c("mat", "desc")]), alpha = fun_alpha,
               permanova = fun_perm, ancombc2 = fun_ct),
          file.path(root, "results.rds"))

  ## Inclusion check: every functional group in every analysis step, condition
  ## and comparison, verified from the files actually written
  incl <- expand_grid(Comparison_id = names(contrasts), Group = fun_groups) |>
    mutate(Condition = map_chr(contrasts[Comparison_id], "cond"),
           Comparison = map_chr(contrasts[Comparison_id], "label"), .before = 1) |>
    mutate(Filtering = Group %in% names(fl),
           Alpha_diversity = file.exists(file.path(d03, paste0("alpha_diversity_", safe_name(Group),
                                                               "_", Condition, ".png"))),
           Beta_diversity_PCoA = file.exists(file.path(d04, paste0("pcoa_bray_", safe_name(Group),
                                                                   "_", Condition, ".png"))),
           Beta_diversity_PCA = file.exists(file.path(d04, paste0("pca_", safe_name(Group),
                                                                  "_", Condition, ".png"))),
           PERMANOVA = paste(Condition, Group) %in%
             paste(fun_perm$Condition, fun_perm$Group)[fun_perm$Comparison == "All groups"],
           ANCOMBC2 = map2_lgl(Comparison_id, Group, ~ .y %in% fun_ct[[.x]]$all$Group),
           Top_features = file.exists(file.path(d06, Comparison_id, paste0("top", cfg$top_n_da, "_features"),
                                                paste0("top", cfg$top_n_da, "_", safe_name(Group), ".png"))))
  save_table(incl, root, "functional_groups_included_by_step")
  incomplete <- incl |> dplyr::filter(!if_all(Filtering:Top_features))
  if (nrow(incomplete) > 0) {
    log_msg("WARNING: functional group x comparison combinations missing from at least one step: ",
            paste(incomplete$Group, incomplete$Comparison_id, collapse = "; "),
            " (see functional_groups_included_by_step.csv)")
  } else {
    log_msg("All ", length(fun_groups), " functional groups analysed in every step for every condition and comparison")
  }
  log_msg("Functional analysis written to ", root)
}

## ---- 7. Run ------------------------------------------------------------------

results_taxa <- try_step("Taxonomic analysis", analyse_taxa())
try_step("Functional analysis", analyse_functions())

writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
log_msg("Done. Outputs in ", normalizePath(cfg$out_dir))
