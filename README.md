# Shotgun metagenome analysis in R (nf-core/mag → phyloseq)

A tutorial for analysing shotgun metagenomes processed with nf-core/mag 5.5.0
and eggNOG-mapper, in R with phyloseq, microViz, microbiome, vegan and ANCOMBC.
The example data are 30 canine faecal metagenomes with three conditions:
gender (female, male), fecal_color (Normal, Abnormally Dark) and fecal_odor
(Normal, Moderately Unpleasant, Unpleasant).

**[Open the tutorial](01_shotgun_tutorial.md)**

| Section | Content |
|---|---|
| 1 | Upstream processing on the MSU HPCC: nf-core/mag (QC, host removal, assembly, binning, DAS Tool, CheckM2, GTDB-Tk), functional annotation with eggNOG-mapper, phyloseq objects |
| 2–11 | GTDB-Tk MAG taxonomy: design, raw abundance, library size, rarefaction, alpha diversity, relative abundance, PCoA, PERMANOVA, ANCOM-BC2 with top-5 boxplots, summary |
| 12 | Functional analysis of COG, EC, Pfam, KEGG KO, KEGG Module, KEGG Pathway: alpha diversity, PCoA, PERMANOVA, ANCOM-BC2 with top-5 boxplots, summary |

## Repository layout

```
01_shotgun_tutorial.Rmd   source (R Markdown)
01_shotgun_tutorial.md    rendered page shown on GitHub
data/                     sample metadata
scripts/                  HPCC and R scripts (Section 1)
results/                  result tables
figures/                  figures (gtdbtk/, function/)
```

## Scripts

| Script | Step |
|---|---|
| `scripts/prep_references.sh` | dog + human host reference (once) |
| `scripts/make_samplesheet.py` | nf-core/mag samplesheet |
| `scripts/run_mag.slurm` | nf-core/mag 5.5.0 SLURM job |
| `scripts/func_config.sh`, `scripts/setup_eggnog_taxa.sh` | functional annotation settings and eggNOG taxonomy files |
| `scripts/1_prepare.sh` … `scripts/6_name_features.sh` | functional annotation, bacteria-only filter, names |
| `scripts/7_create_phyloseq.R` | phyloseq objects |
| `scripts/8_metagenome_analysis.R` | taxonomic and functional analysis |

## Re-render

```r
install.packages(c("rmarkdown", "knitr", "readr", "dplyr", "tidyr"))
rmarkdown::render("01_shotgun_tutorial.Rmd", output_format = "github_document")
```
