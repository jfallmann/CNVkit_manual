required_pkgs <- c("dplyr","ggplot2","purrr","readr","stringr","tidyr","jsonlite",
                   "scales","ggrepel","ggalluvial","patchwork","plotly","htmlwidgets")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Required R packages missing: ", paste(missing_pkgs, collapse=", "))

library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(jsonlite)

manifest_path <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(manifest_path) || !file.exists(manifest_path))
  stop("Usage: Rscript cnv_summary.R <manifest.json> (not found: ", manifest_path, ")")
manifest <- fromJSON(manifest_path, simplifyVector = TRUE)

OUT_DIR <- manifest$out_dir
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TH <- manifest$thresholds
LOG2_HIGH   <- TH$log2_high
LOG2_LOW    <- TH$log2_low
PROBE_HIGH  <- TH$probe_high
PROBE_MED   <- TH$probe_med
MIN_SIZE_KB <- TH$min_size_kb
GAP_MERGE   <- TH$gap_merge
GENE_LOG2   <- TH$gene_log2
GENE_WEIGHT <- TH$gene_weight

PL <- manifest$plots
TOP_N       <- PL$top_n
HM_MAX      <- PL$heatmap_max
CN_CAP      <- PL$cn_cap
HEATMAP_CAP <- PL$heatmap_cap
LABEL_TOP   <- PL$label_top
LOG2_CAP    <- PL$log2_cap

ONCOGENES <- manifest$oncogenes
RESISTANCE_SAMPLES <- manifest$resistance_samples

comparisons <- manifest$comparisons
comparison_ids <- comparisons$id
comparison_type_of <- setNames(comparisons$comparison_type, comparisons$id)
ploidy_of <- setNames(comparisons$ploidy, comparisons$id)
purity_of <- setNames(comparisons$purity, comparisons$id)
sample_of <- setNames(comparisons$sample, comparisons$id)
reference_of <- setNames(comparisons$reference, comparisons$id)
mode_of <- setNames(comparisons$mode, comparisons$id)
cnr_of <- setNames(comparisons$cnr, comparisons$id)
cns_of <- setNames(comparisons$cns, comparisons$id)
gm_of <- setNames(comparisons$genemetrics, comparisons$id)

expected_outputs <- manifest$expected_outputs

GENEMETRICS_FILES <- list()
for (i in seq_len(nrow(comparisons))) {
  s <- comparisons$sample[i]
  if (is.null(GENEMETRICS_FILES[[s]])) GENEMETRICS_FILES[[s]] <- comparisons$genemetrics[i]
}

TEMPORAL_CLONES <- list()
if (length(manifest$temporal_clones) > 0) {
  TEMPORAL_CLONES <- lapply(seq_len(nrow(manifest$temporal_clones)), function(i) {
    spec <- manifest$temporal_clones[i, ]
    list(clone = spec$clone, out_subdir = spec$out_subdir,
         timeline = as.data.frame(spec$timeline[[1]]))
  })
}

CONSERVED_COMPARISONS <- list()
if (length(manifest$conserved_comparisons) > 0) {
  CONSERVED_COMPARISONS <- lapply(seq_len(nrow(manifest$conserved_comparisons)), function(i) {
    spec <- manifest$conserved_comparisons[i, ]
    list(clones = spec$clones[[1]], out_subdir = spec$out_subdir)
  })
}

SCATTER_COMPARISONS <- list()
if (length(manifest$scatter_comparisons) > 0) {
  SCATTER_COMPARISONS <- lapply(seq_len(nrow(manifest$scatter_comparisons)), function(i) {
    spec <- manifest$scatter_comparisons[i, ]
    sc <- list(out_prefix = spec$out_prefix, title = spec$title,
               xlab = spec$xlab, ylab = spec$ylab,
               timeline = as.data.frame(spec$timeline[[1]]))
    if (!is.null(spec$gm) && !is.na(spec$gm$x)) sc$gm <- list(x = spec$gm$x, y = spec$gm$y)
    sc
  })
}



# -- 2. Helper: parse one comparison -------------------------------------------
purity_val <- function(cid) {
  p <- purity_of[[cid]]
  if (is.null(p)) NA_real_ else as.numeric(p)
}

# Read a typed TSV with required-column validation.
read_typed_tsv <- function(path, required, numeric_cols) {
  if (!file.exists(path)) stop("missing input file: ", path)
  df <- read_tsv(path, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    rename_with(tolower)
  missing <- setdiff(required, names(df))
  if (length(missing) > 0)
    stop("missing required column(s) in ", path, ": ", paste(missing, collapse=", "))
  for (col in intersect(numeric_cols, names(df))) {
    vals <- df[[col]]
    num <- suppressWarnings(as.numeric(vals))
    bad <- is.na(num) & !is.na(vals) & !(trimws(vals) %in% c("", "NA", "na", "NaN", "nan"))
    if (any(bad)) {
      i <- which(bad)[1]
      stop("malformed numeric value in column ", col, " of ", path, ": ", vals[i])
    }
    df[[col]] <- num
  }
  df
}

# Extract clean HGNC-style gene symbols from a .cns concatenated gene string
# (segments list transcripts/IDs separated by ';' and genes by ',').
extract_symbols <- function(gene_str) {
  toks <- unlist(str_split(gene_str, "[,;]"))
  toks <- toks[str_detect(toks, "^[A-Z][A-Z0-9-]+$")]
  toks <- toks[!str_detect(toks,
    "^(ENSG|ENST|ENSP|NM|NR|XM|XR|NP|YP|LOC|CCDS|LINC|MIR|SNOR|RNA|ClinID|AC[0-9]|AL[0-9]|AP[0-9]|AF[0-9]|Z[0-9]|U[0-9])")]
  unique(toks)
}

parse_comparison <- function(cid) {
  comparison <- cid
  meta <- list(comparison_type = comparison_type_of[[cid]],
               tumor = sample_of[[cid]], reference = reference_of[[cid]])

  # -- 2a. Read CNR (bin-level log2 ratios, for background scatter) ----------
  cnr_file <- cnr_of[[cid]]
  if (!file.exists(cnr_file)) stop("missing .cnr for ", cid, ": ", cnr_file)
  cnr <- read_typed_tsv(cnr_file,
                        required = c("chromosome","start","end","log2"),
                        numeric_cols = c("start","end","log2")) %>%
    mutate(comparison = comparison)

  # -- 2b. Read called segments (CNS) ----------------------------------------
  cns_file <- cns_of[[cid]]
  if (!file.exists(cns_file)) stop("missing .cns for ", cid, ": ", cns_file)

  cns_raw <- read_typed_tsv(cns_file,
                            required = c("chromosome","start","end","gene","log2","cn","probes"),
                            numeric_cols = c("start","end","log2","cn","probes","weight","depth"))

  seg <- cns_raw %>%
    mutate(
      chromosome      = as.character(chromosome),
      start           = as.numeric(start),
      end             = as.numeric(end),
      svlen           = end - start,
      probes          = if ("probes" %in% names(.)) as.numeric(probes) else NA_real_,
      cn              = if ("cn" %in% names(.)) as.numeric(cn) else NA_real_,
      log2fc          = as.numeric(log2),
      # the .cns gene column concatenates every gene in the segment and can be
      # tens of thousands of characters; keep only a count to keep tables lean
      n_genes         = if ("gene" %in% names(.)) str_count(gene, ",") + 1L else NA_integer_,
      comparison      = comparison,
      comparison_type = meta$comparison_type,
      tumor           = meta$tumor,
      reference       = meta$reference,
      id              = cid,
      mode            = mode_of[[cid]],
      sample          = sample_of[[cid]],
      ploidy          = ploidy_of[[cid]],
      purity          = purity_val(cid)
    ) %>%
    select(-any_of("gene"))

  # -- 2c. Read genemetrics --------------------------------------------------
  gm_file <- gm_of[[cid]]
  if (!file.exists(gm_file)) stop("missing genemetrics for ", cid, ": ", gm_file)
  gm <- read_typed_tsv(gm_file,
                       required = c("gene","chromosome","start","end","log2","weight","cn"),
                       numeric_cols = c("start","end","probes","log2","weight","cn")) %>%
    mutate(comparison      = comparison,
           comparison_type = meta$comparison_type,
           tumor           = meta$tumor,
           reference       = meta$reference,
           id              = cid,
           mode            = mode_of[[cid]],
           sample          = sample_of[[cid]],
           ploidy          = ploidy_of[[cid]],
           purity          = purity_val(cid))

  list(cnr=cnr, seg=seg, gm=gm)
}

# -- 3. Parse all comparisons --------------------------------------------------
all_data <- map(comparison_ids, parse_comparison)
names(all_data) <- comparison_ids

# -- 4. Trustworthiness scoring ------------------------------------------------
# Score each called segment 0-100 based on:
#   - probe support     (0-40 pts)
#   - |log2FC| magnitude (0-40 pts)
#   - segment size      (0-20 pts)

score_segment <- function(probes, log2fc, svlen_kb) {
  # Probe score (40 pts max)
  probe_score <- case_when(
    probes >= PROBE_HIGH ~ 40,
    probes >= PROBE_MED  ~ 20,
    probes >= 10         ~ 10,
    TRUE                 ~ 0
  )
  # log2FC score (40 pts max)
  fc_score <- case_when(
    abs(log2fc) >= 1.0  ~ 40,
    abs(log2fc) >= 0.6  ~ 30,
    abs(log2fc) >= LOG2_HIGH ~ 20,
    abs(log2fc) >= LOG2_LOW  ~ 10,
    TRUE                ~ 0
  )
  # Size score (20 pts max)
  size_score <- case_when(
    svlen_kb >= 5000 ~ 20,
    svlen_kb >= 1000 ~ 15,
    svlen_kb >= 500  ~ 10,
    svlen_kb >= MIN_SIZE_KB ~ 5,
    TRUE             ~ 0
  )
  probe_score + fc_score + size_score
}

label_confidence <- function(score) {
  case_when(
    score >= 70 ~ "HIGH",
    score >= 40 ~ "MODERATE",
    score >= 20 ~ "LOW",
    TRUE        ~ "NOISE"
  )
}

# Combine all called segments and score
all_seg <- map_dfr(all_data, "seg") %>%
  filter(!is.na(log2fc)) %>%
  mutate(
    svlen_kb   = svlen / 1000,
    probes_eff = ifelse(is.na(probes), PROBE_MED, probes),
    score      = score_segment(probes_eff, log2fc, svlen_kb),
    confidence = label_confidence(score),
    direction  = ifelse(log2fc > 0, "GAIN", "LOSS"),
    # clean chromosome
    chromosome = ifelse(startsWith(chromosome,"chr"), chromosome,
                        paste0("chr", chromosome))
  ) %>%
  filter(abs(log2fc) >= LOG2_LOW, svlen_kb >= MIN_SIZE_KB) %>%
  arrange(desc(score), desc(abs(log2fc)))

# -- 5. Gene-level summary from genemetrics -----------------------------------
all_gm <- map_dfr(all_data, "gm") %>%
  filter(!is.na(log2)) %>%
  mutate(
    direction  = ifelse(log2 > 0, "GAIN", "LOSS"),
    confidence = label_confidence(
      score_segment(
        probes = ifelse("probes" %in% names(.), probes, 50),  # default if missing
        log2fc = log2,
        svlen_kb = 500  # gene-level, size not applicable
      )
    )
  ) %>%
  filter(abs(log2) >= LOG2_LOW) %>%
  arrange(desc(abs(log2)))

# -- 6. Write output tables ----------------------------------------------------
# Split by comparison type: differences vs the human reference (ancestral +
# acquired) vs divergence between populations (pairwise).
seg_ref  <- all_seg %>% filter(comparison_type == "vs_reference")
seg_pair <- all_seg %>% filter(comparison_type == "pairwise")

write_csv(all_seg,  file.path(OUT_DIR, "all_cnv_segments_scored.csv"))
write_csv(seg_ref,  file.path(OUT_DIR, "segments_scored_vs_reference.csv"))
write_csv(seg_pair, file.path(OUT_DIR, "segments_scored_pairwise.csv"))
write_csv(all_seg %>% filter(confidence %in% c("HIGH","MODERATE")),
          file.path(OUT_DIR, "cnv_segments_high_moderate.csv"))
write_csv(all_gm,   file.path(OUT_DIR, "all_genemetrics_combined.csv"))

# One scored-segments file per comparison
by_dir <- file.path(OUT_DIR, "by_comparison")
dir.create(by_dir, showWarnings = FALSE)
walk(comparison_ids, ~ write_csv(all_seg %>% filter(comparison == .x),
                                 file.path(by_dir, paste0(.x, ".csv"))))

message("Segments scored: ", nrow(all_seg),
        " (vs_reference: ", nrow(seg_ref), ", pairwise: ", nrow(seg_pair), ")")
message("HIGH confidence: ", sum(all_seg$confidence == "HIGH"))
message("MODERATE:        ", sum(all_seg$confidence == "MODERATE"))

# -- 7. Manhattan-style plot ---------------------------------------------------
# Compute cumulative x-axis positions per chromosome
chrom_order <- paste0("chr", c(1:22,"X","Y"))

# Use CNR bins for the background scatter (all bins, all comparisons)
all_cnr <- map_dfr(all_data, "cnr") %>%
  filter(chromosome %in% chrom_order |
         paste0("chr",chromosome) %in% chrom_order) %>%
  mutate(chromosome = ifelse(startsWith(chromosome,"chr"), chromosome,
                             paste0("chr", chromosome)),
         chromosome = factor(chromosome, levels=chrom_order))

# Chromosome sizes from data (raw cnr bins, else raw segments, else placeholder)
if (nrow(all_cnr) > 0) {
  chrom_sizes <- all_cnr %>%
    group_by(chromosome) %>%
    summarise(max_pos = max(end, na.rm=TRUE), .groups="drop") %>%
    arrange(match(chromosome, chrom_order)) %>%
    mutate(offset = cumsum(lag(max_pos, default=0)),
           mid    = offset + max_pos/2)
} else if (nrow(all_seg) > 0) {
  chrom_sizes <- all_seg %>%
    group_by(chromosome) %>%
    summarise(max_pos = max(end, na.rm=TRUE), .groups="drop") %>%
    arrange(match(chromosome, chrom_order)) %>%
    mutate(offset = cumsum(lag(max_pos, default=0)),
           mid    = offset + max_pos/2)
} else {
  chrom_sizes <- tibble(chromosome = factor(character(), levels = chrom_order),
                        max_pos = numeric(), offset = numeric(), mid = numeric())
}

# Add cumulative position
all_cnr <- all_cnr %>%
  left_join(chrom_sizes[,c("chromosome","offset")], by="chromosome") %>%
  mutate(x_pos = (start + end)/2 + offset)

all_seg_plot <- all_seg %>%
  mutate(chromosome = factor(chromosome, levels=chrom_order)) %>%
  left_join(chrom_sizes[,c("chromosome","offset")], by="chromosome") %>%
  mutate(
    x_start = start + offset,
    x_mid   = (start + end)/2 + offset,
    x_end   = end + offset,
    # order facets so vs_reference panels precede pairwise panels
    comparison = factor(comparison,
                        levels = unique(comparison[order(comparison_type, comparison)]))
  )

# Alternating chromosome colors for background
chrom_colors <- setNames(
  rep(c("grey80","grey92"), length(chrom_order)),
  chrom_order
)

# One plot per comparison, faceted
p <- ggplot() +
  # Background chromosome bands
  geom_rect(data = chrom_sizes,
            aes(xmin=offset, xmax=offset+max_pos,
                ymin=-Inf, ymax=Inf,
                fill=chromosome),
            alpha=0.3, show.legend=FALSE) +
  scale_fill_manual(values=chrom_colors) +

  # CNR bins - grey background scatter
  geom_point(data = all_cnr,
             aes(x=x_pos, y=log2),
             size=0.08, alpha=0.15, color="grey50") +

  # Called segments as thick horizontal lines, colored by confidence
  geom_segment(data = all_seg_plot,
               aes(x=x_start, xend=x_end,
                   y=log2fc, yend=log2fc,
                   color=confidence),
               linewidth=1.2, alpha=0.9) +

  # Confidence color scale
  scale_color_manual(
    values = c("HIGH"="firebrick","MODERATE"="darkorange",
               "LOW"="steelblue","NOISE"="grey70"),
    name   = "Confidence"
  ) +

  # Reference lines
  geom_hline(yintercept=0,         color="black",  linewidth=0.6) +
  geom_hline(yintercept=c(0.4,-0.4), color="red",  linewidth=0.4, linetype="dashed") +
  geom_hline(yintercept=c(0.2,-0.2), color="blue", linewidth=0.3, linetype="dotted") +

  # Chromosome labels on x-axis
  scale_x_continuous(
    breaks = chrom_sizes$mid,
    labels = str_remove(chrom_sizes$chromosome,"chr"),
    expand = c(0.01, 0)
  ) +

  # Facet by comparison
  {if (nrow(all_seg_plot) > 0)
     facet_wrap(~comparison, ncol=1, strip.position="right") else NULL} +

  {if (nrow(all_seg_plot) == 0)
     annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.1,
              label="No called segments passed filters", size=4) else NULL} +

  # Clip the y-axis to a fixed +/- LOG2_CAP window (shared across all facets,
  # since facet_wrap uses a common scale). A hard cap - rather than a
  # quantile of the data - keeps the axis readable even when a handful of
  # segments (e.g. chrY homozygous deletions at log2fc ~ -7) would otherwise
  # still dominate the 2nd/98th percentile and blow out every panel.
  # Outlier segments are still drawn but clipped at the panel edge.
  coord_cartesian(ylim = c(-LOG2_CAP, LOG2_CAP)) +

  labs(
    x     = "Chromosome",
    y     = "log2 Fold Change (vs reference / comparison)",
    title = "CNV landscape across comparisons",
    subtitle = "Segments colored by confidence | dashed=0.4 | dotted=0.2 threshold"
  ) +
  theme_bw(base_size=11) +
  theme(
    axis.text.x     = element_text(size=7, angle=0),
    panel.grid      = element_blank(),
    strip.text.y    = element_text(size=8, angle=0),
    legend.position = "bottom"
  )

ggsave(file.path(OUT_DIR, "manhattan_cnv_all_comparisons.pdf"),
       p, width=18, height=4*length(comparison_ids), limitsize=FALSE)
ggsave(file.path(OUT_DIR, "manhattan_cnv_all_comparisons.png"),
       p, width=18, height=4*length(comparison_ids), dpi=150, limitsize=FALSE)

message("Plot saved to ", OUT_DIR)

# -- 8. Comparison-wise region matrices ---------------------------------------
# Collapse overlapping, same-direction segments into regions, spread the
# per-comparison log2FC across columns, and classify each region as
# comparison-specific or conserved across the comparisons in that set.
CONF_LEVELS <- c("HIGH","MODERATE","LOW","NOISE")

build_region_matrix <- function(df, all_ids) {
  n_total <- length(all_ids)
  if (nrow(df) == 0) {
    out <- tibble(chromosome = character(), direction = character(), region_id = integer(),
                  region_start = integer(), region_end = integer(), n_comparisons = integer(),
                  comparisons = character(), mean_log2fc = numeric(), max_score = numeric(),
                  max_confidence = character(), class = character())
    for (cid in all_ids) out[[cid]] <- numeric()
    return(out)
  }
  regions <- df %>%
    arrange(chromosome, direction, start) %>%
    group_by(chromosome, direction) %>%
    mutate(region_id = cumsum(start > lag(cummax(end), default = 0) + GAP_MERGE)) %>%
    ungroup()

  region_meta <- regions %>%
    group_by(chromosome, direction, region_id) %>%
    summarise(
      region_start   = min(start),
      region_end     = max(end),
      n_comparisons  = n_distinct(comparison),
      comparisons    = paste(sort(unique(comparison)), collapse="; "),
      mean_log2fc    = mean(log2fc),
      max_score      = max(score),
      max_confidence = CONF_LEVELS[min(match(confidence, CONF_LEVELS))],
      .groups="drop"
    ) %>%
    mutate(
      class = case_when(
        n_total > 1 & n_comparisons == n_total ~ "conserved_all",
        n_comparisons == 1                     ~ "comparison_specific",
        TRUE                                   ~ "shared_subset"
      )
    )

  wide <- regions %>%
    group_by(chromosome, direction, region_id, comparison) %>%
    summarise(log2fc = mean(log2fc), .groups="drop") %>%
    pivot_wider(names_from = comparison, values_from = log2fc, values_fill = 0)
  for (cid in all_ids) if (!cid %in% names(wide)) wide[[cid]] <- 0

  region_meta %>%
    left_join(wide, by = c("chromosome","direction","region_id")) %>%
    arrange(desc(max_score))
}

# Does a region overlap a same-direction reference region?
overlaps_reference <- function(chrom, dir, s, e, ref) {
  if (is.null(ref) || nrow(ref) == 0) return(FALSE)
  any(ref$chromosome == chrom & ref$direction == dir &
      ref$region_start <= e & ref$region_end >= s)
}

hm <- function(df) df %>% filter(confidence %in% c("HIGH","MODERATE"))

ref_ids  <- comparison_ids[comparison_type_of[comparison_ids] == "vs_reference"]
pair_ids <- comparison_ids[comparison_type_of[comparison_ids] == "pairwise"]

regions_ref  <- build_region_matrix(hm(seg_ref), ref_ids)
regions_pair <- build_region_matrix(hm(seg_pair), pair_ids)

# Annotate pairwise regions with whether they overlap a same-direction
# vs_reference region (overlap heuristic only, not evidence of ancestry).
regions_pair <- regions_pair %>%
  rowwise() %>%
  mutate(overlaps_vs_reference = overlaps_reference(chromosome, direction,
                                 region_start, region_end, regions_ref)) %>%
  ungroup()

write_csv(regions_ref, file.path(OUT_DIR, "cnv_regions_vs_reference.csv"))
message("vs_reference regions: ", nrow(regions_ref),
        " (conserved_all: ", sum(regions_ref$class == "conserved_all"),
        ", comparison_specific: ", sum(regions_ref$class == "comparison_specific"), ")")

write_csv(regions_pair, file.path(OUT_DIR, "cnv_regions_pairwise.csv"))
message("pairwise regions: ", nrow(regions_pair),
        " (conserved_all: ", sum(regions_pair$class == "conserved_all"),
        ", comparison_specific: ", sum(regions_pair$class == "comparison_specific"),
        "; overlapping vs_reference: ", sum(regions_pair$overlaps_vs_reference), ")")

# Slim table of regions conserved across ALL pairwise comparisons.
shared_all <- regions_pair %>%
  filter(class == "conserved_all") %>%
  dplyr::select(chromosome, direction, region_id, region_start, region_end,
                n_comparisons, comparisons, mean_log2fc, max_score)
write_csv(shared_all, file.path(OUT_DIR, "cnv_regions_shared_all_comparisons.csv"))
message("shared-all-comparison regions: ", nrow(shared_all))


# ==============================================================================
# 8. Gene-level genemetrics summary (config-driven; skipped when no files present)
# ==============================================================================
# Builds two wide tables of gene-level log2 CNVs from the validated all_gm
# table, with one column per configured comparison ID (all baselines retained).
# Samples named in RESISTANCE_SAMPLES drive the resistance-specific filter.

gm_base <- all_gm %>%
  filter(abs(log2) > GENE_LOG2, weight > GENE_WEIGHT) %>%
  mutate(sym = map(gene, extract_symbols)) %>%
  tidyr::unnest(sym) %>%
  filter(!is.na(sym), nzchar(sym))

# Pivot to wide - one row per gene, one col per comparison ID
gm_wide <- gm_base %>%
  group_by(sym, chromosome, id) %>%
  summarise(log2 = mean(log2, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = id, values_from = log2, values_fill = NA)
for (cid in comparison_ids) if (!cid %in% names(gm_wide)) gm_wide[[cid]] <- NA_real_
gm_wide <- gm_wide %>%
  dplyr::select(sym, chromosome, all_of(comparison_ids)) %>%
  rename(gene = sym)

# Genes altered in resistant clones (comparison IDs whose sample is resistant)
resist_ids <- comparison_ids[sample_of[comparison_ids] %in% RESISTANCE_SAMPLES]
resistance_specific <- if (length(resist_ids) > 0) {
  gm_wide %>%
    filter(if_any(all_of(resist_ids), ~ abs(.x) > GENE_LOG2))
} else {
  gm_wide[0, ]
}

write.table(gm_wide,
            file.path(OUT_DIR, "gene_log2_matrix.tsv"),
            sep="\t", row.names = F, quote = F)

write.table(resistance_specific,
            file.path(OUT_DIR, "resistant_group_genes.tsv"),
            sep="\t", row.names = F, quote = F)

# Clean gene names and summarise per primary gene symbol
gm_clean <- all_gm %>%
  mutate(sym = map(gene, extract_symbols)) %>%
  tidyr::unnest(sym) %>%
  filter(!is.na(sym), nzchar(sym), weight > GENE_WEIGHT) %>%
  group_by(sym, chromosome, id) %>%
  summarise(
    region_start = min(start, na.rm=TRUE),
    region_end   = max(end,   na.rm=TRUE),
    log2_mean    = mean(log2, na.rm=TRUE),
    log2_max     = max(log2,  na.rm=TRUE),
    log2_min     = min(log2,  na.rm=TRUE),
    log2_extreme = ifelse(abs(max(log2)) > abs(min(log2)), max(log2), min(log2)),
    n_bins       = n(),
    .groups      = "drop"
  )

# Pivot log2_extreme wide; keep coordinates by taking min/max across comparisons
coords <- gm_clean %>%
  group_by(sym, chromosome) %>%
  summarise(
    start = min(region_start, na.rm=TRUE),
    end   = max(region_end,   na.rm=TRUE),
    .groups = "drop"
  )

gm_wide_clean <- gm_clean %>%
  filter(abs(log2_extreme) > GENE_LOG2) %>%
  dplyr::select(sym, chromosome, id, log2_extreme) %>%
  pivot_wider(
    names_from  = id,
    values_from = log2_extreme,
    values_fill = NA
  )
for (cid in comparison_ids) if (!cid %in% names(gm_wide_clean)) gm_wide_clean[[cid]] <- NA_real_
gm_wide_clean <- gm_wide_clean %>%
  left_join(coords, by = c("sym", "chromosome")) %>%
  dplyr::select(sym, chromosome, start, end, all_of(comparison_ids)) %>%
  arrange(chromosome, start) %>%
  rename(gene = sym)

write.table(gm_wide_clean,
            file.path(OUT_DIR, "gene_log2_matrix_cleaned.tsv"),
            sep="\t", row.names = F, quote = F)

# ==============================================================================
# 9. Per-clone temporal copy-number drift (BC139, BC217)
# ==============================================================================
# Runs the same lineage analysis for each clone: a relative log2 Manhattan, a
# gene x timepoint CN heatmap, per-gene CN trajectories and a CN-state alluvial.
# BC139 has a resistant derivative (4th timepoint); BC217 does not (3 timepoints).
# Absolute copy number is taken from the called-segment 'cn' column of the
# .call.cns (preferring somatic.call.cns via find_cns()).

# Read one timepoint's called segments and explode to one row per gene symbol.
read_timepoint <- function(timepoint, label, dir) {
  f <- cns_of[[dir]]
  if (is.null(f) || is.na(f)) { message("No .cns for timepoint ", label, " (", dir, ") - skipped"); return(NULL) }
  if (!file.exists(f)) stop("missing .cns for timepoint ", label, " (", dir, "): ", f)
  read_typed_tsv(f, required = c("chromosome","start","end","gene","log2","cn","probes"),
                 numeric_cols = c("start","end","log2","cn","probes","weight","depth")) %>%
    filter(!is.na(cn)) %>%
    mutate(chromosome = as.character(chromosome),
           sym = map(gene, extract_symbols)) %>%
    tidyr::unnest(sym) %>%
    transmute(timepoint, label, gene = sym, chromosome, start, end, log2, cn,
              ploidy = ploidy_of[[dir]], comparison = dir,
              reference = reference_of[[dir]], purity = purity_val(dir),
              sample = sample_of[[dir]])
}

# Collapse a character vector into a hover-friendly list, capping the length so
# a coordinate shared by thousands of genes does not produce a giant tooltip.
fmt_gene_list <- function(x, k = 60) {
  if (length(x) > k)
    paste0(paste(x[seq_len(k)], collapse = "<br>"), "<br>... +", length(x) - k, " more")
  else paste(x, collapse = "<br>")
}

# Interactive plotly companion to the per-gene CN scatters. Points are aggregated
# by (clamped x, clamped y, state); dot area encodes the number of genes and the
# hover lists every gene at that coordinate (with its chr:start locus when known),
# so overlapping same-CN genes are browsable. Writes <out_prefix>.html. Skips
# gracefully when plotly/htmlwidgets are not installed.
save_scatter_html <- function(df, xlab, ylab, title, out_prefix, state_cols, cap = CN_CAP) {
  if (nrow(df) == 0) {
    p <- plotly::plot_ly() %>%
      plotly::layout(title = title,
                     xaxis = list(title = xlab, range = c(-0.5, cap + 0.5)),
                     yaxis = list(title = ylab, range = c(-0.5, cap + 0.5)),
                     annotations = list(text = "No shared genes", showarrow = FALSE,
                                        x = 0.5, y = 0.5, xref = "paper", yref = "paper"))
    htmlwidgets::saveWidget(p, paste0(out_prefix, ".html"), selfcontained = TRUE)
    message("Interactive scatter written: ", basename(out_prefix), ".html")
    return(invisible(NULL))
  }
  has_loc <- all(c("chr", "start") %in% names(df))
  df <- df %>%
    mutate(xp = pmin(x, cap), yp = pmin(y, cap),
           locus = if (has_loc) paste0(gene, " (", chr, ":", start, ")") else gene)
  agg <- df %>%
    group_by(xp, yp, state) %>%
    summarise(n = n(), genes = fmt_gene_list(locus), .groups = "drop") %>%
    mutate(hover = paste0("CN: (", xp, ", ", yp, ")<br>state: ", state,
                          "<br>", n, " gene(s)<br>", genes))
  p <- plotly::plot_ly(agg, x = ~xp, y = ~yp, type = "scatter", mode = "markers",
                       color = ~state, colors = state_cols, size = ~n,
                       sizes = c(20, 500),
                       marker = list(sizemode = "area", opacity = 0.6),
                       text = ~hover, hoverinfo = "text") %>%
    plotly::layout(title = title,
                   xaxis = list(title = xlab, range = c(-0.5, cap + 0.5)),
                   yaxis = list(title = ylab, range = c(-0.5, cap + 0.5)))
  htmlwidgets::saveWidget(p, paste0(out_prefix, ".html"), selfcontained = TRUE)
  message("Interactive scatter written: ", basename(out_prefix), ".html")
}

# Per-gene CN scatter of two comparisons (x vs y). Gains/losses are called
# relative to each axis's own modal CN so a global ploidy shift is not mistaken
# for gene-level change; genes moving the same way on both axes are "in both".
# 'tab' needs columns gene, x, y, is_oncogene (chromosome/start optional, used for
# the interactive hover). Writes <out_prefix>.{pdf,png,html}.
cn_state_scatter <- function(tab, xlab, ylab, title, out_prefix, cap = CN_CAP,
                             label_top = LABEL_TOP, bx = NULL, by = NULL) {
  if (nrow(tab) == 0) {
    p <- ggplot() +
      annotate("text", x = cap/2, y = cap/2, label = "No shared genes", size = 4) +
      coord_cartesian(xlim = c(0, cap), ylim = c(0, cap)) +
      labs(x = xlab, y = ylab, title = title) +
      theme_bw(base_size = 10)
    ggsave(paste0(out_prefix, ".pdf"), p, width = 7, height = 6, limitsize = FALSE)
    ggsave(paste0(out_prefix, ".png"), p, width = 7, height = 6, dpi = 150, limitsize = FALSE)
    save_scatter_html(tab, xlab, ylab, title, out_prefix,
                      c("gain in both" = "#B2182B"), cap = cap)
    return(invisible(tab))
  }
  if (is.null(bx)) bx <- round(median(tab$x, na.rm = TRUE))
  if (is.null(by)) by <- round(median(tab$y, na.rm = TRUE))
  st <- function(cn, base) case_when(cn >= base + 1 ~ "gain",
                                     cn <= base - 1 ~ "loss",
                                     TRUE           ~ "neutral")
  tab <- tab %>%
    mutate(state_x = st(x, bx), state_y = st(y, by),
           joint = case_when(
             state_x == "gain"    & state_y == "gain"    ~ "gain in both",
             state_x == "loss"    & state_y == "loss"    ~ "loss in both",
             state_x == "neutral" & state_y == "neutral" ~ "neutral in both",
             TRUE                                         ~ "divergent"),
           xp = pmin(x, cap), yp = pmin(y, cap))
  # Concordant genes eligible for labelling. Oncogenes are always labelled;
  # the strongest-changing non-oncogenes (top `label_top` by combined distance
  # from the per-axis modal CN) are also labelled so notable genes outside the
  # curated list are not missed.
  concordant <- tab %>% filter(joint %in% c("gain in both", "loss in both"))
  lab_onco  <- concordant %>% filter(is_oncogene)
  lab_other <- concordant %>%
    filter(!is_oncogene) %>%
    mutate(dev = abs(x - bx) + abs(y - by)) %>%
    slice_max(dev, n = label_top, with_ties = FALSE)
  # Divergent genes (changed on one axis but not the other, or opposite
  # direction) are often the most interesting - they are what differs
  # between the two comparisons. Oncogenes are always labelled; the
  # strongest-diverging non-oncogenes (top `label_top` by how far apart the
  # two axes' per-gene changes are) are also labelled.
  divergent <- tab %>% filter(joint == "divergent")
  lab_onco_div  <- divergent %>% filter(is_oncogene)
  lab_other_div <- divergent %>%
    filter(!is_oncogene) %>%
    mutate(dev = abs((x - bx) - (y - by))) %>%
    slice_max(dev, n = label_top, with_ties = FALSE)
  state_cols <- c("gain in both" = "#B2182B", "loss in both" = "#2166AC",
                  "divergent" = "grey60", "neutral in both" = "grey85")
  p <- ggplot(tab, aes(xp, yp, color = joint)) +
    geom_count(alpha = 0.6) +
    scale_size_area(max_size = 9, name = "genes") +
    geom_vline(xintercept = bx, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = by, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = state_cols, name = "State") +
    {if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab_other) > 0)
       ggrepel::geom_text_repel(
         data = lab_other,
         aes(label = gene), size = 2.2, color = "grey25", fontface = "plain",
         max.overlaps = Inf, show.legend = FALSE, inherit.aes = TRUE,
         min.segment.length = 0, segment.size = 0.2, segment.color = "grey70",
         box.padding = 0.5, point.padding = 0.3, force = 3, seed = 1,
         bg.color = "white", bg.r = 0.15)} +
    {if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab_onco) > 0)
       ggrepel::geom_text_repel(
         data = lab_onco,
         aes(label = gene), size = 2.6, color = "black", fontface = "bold",
         max.overlaps = Inf, show.legend = FALSE, inherit.aes = TRUE,
         min.segment.length = 0, segment.size = 0.3, segment.color = "grey40",
         box.padding = 0.6, point.padding = 0.3, force = 3, seed = 1,
         bg.color = "white", bg.r = 0.18)} +
    {if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab_other_div) > 0)
       ggrepel::geom_text_repel(
         data = lab_other_div,
         aes(label = gene), size = 2.2, color = "grey40", fontface = "italic",
         max.overlaps = Inf, show.legend = FALSE, inherit.aes = TRUE,
         min.segment.length = 0, segment.size = 0.2, segment.color = "grey70",
         box.padding = 0.5, point.padding = 0.3, force = 3, seed = 1,
         bg.color = "white", bg.r = 0.15)} +
    {if (requireNamespace("ggrepel", quietly = TRUE) && nrow(lab_onco_div) > 0)
       ggrepel::geom_text_repel(
         data = lab_onco_div,
         aes(label = gene), size = 2.6, color = "black", fontface = "bold.italic",
         max.overlaps = Inf, show.legend = FALSE, inherit.aes = TRUE,
         min.segment.length = 0, segment.size = 0.3, segment.color = "grey40",
         box.padding = 0.6, point.padding = 0.3, force = 3, seed = 1,
         bg.color = "white", bg.r = 0.18)} +
    coord_cartesian(xlim = c(0, cap), ylim = c(0, cap)) +
    labs(x = xlab, y = ylab, title = title,
         subtitle = paste0("Dashed = configured ploidy (", bx, ", ", by, ") | axes clamped at ", cap,
                           " | concordant oncogenes (bold) + top ", label_top,
                           " other concordant genes labelled",
                           " | divergent oncogenes (bold italic) + top ", label_top,
                           " other divergent genes labelled (italic)")) +
    theme_bw(base_size = 10)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 7, height = 6, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), p, width = 7, height = 6, dpi = 150, limitsize = FALSE)
  save_scatter_html(tab %>% mutate(state = joint), xlab, ylab, title, out_prefix,
                    state_cols, cap = cap)
  invisible(tab)
}

# Full temporal analysis for one clone lineage. 'timeline' is a data.frame with
# columns timepoint,label,dir; outputs are written under OUT_DIR/out_subdir.
# Returns per-gene CN summaries used by the conserved-CNV comparison (Section 10).
run_clone_temporal <- function(clone, timeline, out_subdir) {
  TEMP_DIR <- file.path(OUT_DIR, out_subdir)
  dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  labels_ord <- timeline$label
  slug <- tolower(clone)
  base_ploidy <- ploidy_of[[timeline$dir[1]]]
  label_ploidy <- setNames(ploidy_of[timeline$dir], labels_ord)
  label_purity <- setNames(vapply(timeline$dir, purity_val, numeric(1)), labels_ord)

# -- 9a. Relative log2 Manhattan across the clone timeline --------------------
# Absolute cn mixes quantities across comparisons: the vs-reference call is true
# integer CN (near-triploid, baseline 3) while the pairwise somatic calls are
# already relative (log2 vs their comparator, baseline 0). Plotting log2fc puts
# every comparison on the SAME neutral-at-0 scale, so gains/losses are directly
# comparable without tracking each ploidy. y is symmetric and clipped so a few
# deep deletions don't squash the rest; strong changes are emphasised.
LOG2_GAIN_HL <- LOG2_HIGH   # emphasise (thicken) segments beyond this |log2|

cn_pal <- setNames(scales::hue_pal()(nrow(timeline)), labels_ord)

cn_seg <- all_seg %>%
  filter(comparison %in% timeline$dir, !is.na(log2fc)) %>%
  left_join(timeline[, c("dir","label")], by = c("comparison" = "dir")) %>%
  mutate(chromosome = factor(chromosome, levels = chrom_order)) %>%
  filter(!is.na(chromosome)) %>%
  left_join(chrom_sizes[, c("chromosome","offset")], by = "chromosome") %>%
  mutate(
    x_start   = start + offset,
    x_end     = end   + offset,
    label     = factor(label, levels = labels_ord),
    highlight = abs(log2fc) >= LOG2_GAIN_HL
  )

p_cn_abs <- NULL
if (nrow(cn_seg) == 0) {
  message(clone, " log2-Manhattan skipped: none of the timeline dirs are in all_seg.")
  p_cn <- ggplot() +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = "No called segments passed filters", size = 4) +
    labs(x = "Chromosome",
         y = expression(log[2]~"ratio (relative to comparator)"),
         title = paste0(clone, " relative copy-number landscape across the ordered comparison series")) +
    theme_bw(base_size = 11)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_log2_relative.pdf")),
         p_cn, width = 18, height = 6, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_log2_relative.png")),
         p_cn, width = 18, height = 6, dpi = 150, limitsize = FALSE)
  p_cn_abs <- ggplot() +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
             label = "No called segments passed filters", size = 4) +
    labs(x = "Chromosome", y = "Called copy number (cn)",
         title = paste0(clone, " called copy-number landscape across the ordered comparison series")) +
    theme_bw(base_size = 11)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_cn_absolute.pdf")),
         p_cn_abs, width = 18, height = 6, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_cn_absolute.png")),
         p_cn_abs, width = 18, height = 6, dpi = 150, limitsize = FALSE)
} else {
  # CNVKit floors homozygous deletions (cn 0) anywhere down to log2 ~ -27, a
  # flooring artifact rather than real gradation. Size the symmetric window from
  # the gain side (+ a threshold floor); deep losses clamp at the bottom edge.
  cn_cap <- max(quantile(cn_seg$log2fc[cn_seg$log2fc > 0], 0.995, na.rm = TRUE),
                LOG2_GAIN_HL * 4, na.rm = TRUE)

  p_cn <- ggplot() +
    geom_rect(data = chrom_sizes,
              aes(xmin = offset, xmax = offset + max_pos, ymin = -Inf, ymax = Inf,
                  fill = chromosome),
              alpha = 0.3, show.legend = FALSE) +
    scale_fill_manual(values = chrom_colors) +
    # shade the gain / loss bands beyond the highlight threshold
    annotate("rect", xmin = -Inf, xmax = Inf, ymin =  LOG2_GAIN_HL, ymax =  Inf,
             fill = "firebrick", alpha = 0.05) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -LOG2_GAIN_HL,
             fill = "steelblue", alpha = 0.05) +
    # references: neutral (log2 0, no change vs comparator) and gain/loss cutoffs
    geom_hline(yintercept =  0,            color = "black",     linewidth = 0.4, linetype = "solid") +
    geom_hline(yintercept =  LOG2_GAIN_HL, color = "firebrick", linewidth = 0.4, linetype = "dashed") +
    geom_hline(yintercept = -LOG2_GAIN_HL, color = "steelblue", linewidth = 0.4, linetype = "dashed") +
    # faint background: segments within the neutral band
    geom_segment(data = dplyr::filter(cn_seg, !highlight),
                 aes(x = x_start, xend = x_end, y = log2fc, yend = log2fc, color = label),
                 linewidth = 0.6, alpha = 0.45) +
    # emphasised foreground: gains/losses beyond the threshold
    geom_segment(data = dplyr::filter(cn_seg, highlight),
                 aes(x = x_start, xend = x_end, y = log2fc, yend = log2fc, color = label),
                 linewidth = 1.8, alpha = 1) +
    scale_color_manual(values = cn_pal, name = "Timepoint") +
    scale_x_continuous(breaks = chrom_sizes$mid,
                       labels = str_remove(chrom_sizes$chromosome, "chr"),
                       expand = c(0.01, 0)) +
    scale_y_continuous(breaks = scales::breaks_width(1)) +
    coord_cartesian(ylim = c(-cn_cap, cn_cap)) +
    labs(
      x = "Chromosome",
      y = expression(log[2]~"ratio (relative to comparator)"),
      title = paste0(clone, " relative copy-number landscape across the ordered comparison series"),
      subtitle = paste0("Single panel, colored by timepoint | solid black = no change (log2 0) | ",
                        "dashed = +/- ", LOG2_GAIN_HL, " gain/loss cutoff | thick = |log2| >= ",
                        LOG2_GAIN_HL, " | y clipped at +/-", round(cn_cap, 2))
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x     = element_text(size = 7),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )

  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_log2_relative.pdf")),
         p_cn, width = 18, height = 6, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_log2_relative.png")),
         p_cn, width = 18, height = 6, dpi = 150, limitsize = FALSE)
  message("log2 Manhattan saved (y clipped at +/-", round(cn_cap, 2), ")")

  # -- absolute copy-number Manhattan (integer CN, baseline = modal ploidy) ----
  # Same layout as the log2 panel but on the raw integer-CN axis, so the
  # near-triploid baseline and focal amplifications are read directly.
  cn_abs_cap <- CN_CAP
  cn_seg_abs <- cn_seg %>% mutate(hl_abs = cn >= (ploidy + 1) | cn <= (ploidy - 1))
  baseline_df <- cn_seg_abs %>%
    distinct(label, ploidy) %>%
    arrange(label)
  p_cn_abs <- ggplot() +
    geom_rect(data = chrom_sizes,
              aes(xmin = offset, xmax = offset + max_pos, ymin = -Inf, ymax = Inf,
                  fill = chromosome),
              alpha = 0.3, show.legend = FALSE) +
    scale_fill_manual(values = chrom_colors) +
    geom_hline(data = baseline_df, aes(yintercept = ploidy, color = label),
               linewidth = 0.4, linetype = "dashed", alpha = 0.8) +
    geom_segment(data = dplyr::filter(cn_seg_abs, !hl_abs),
                 aes(x = x_start, xend = x_end, y = cn, yend = cn, color = label),
                 linewidth = 0.6, alpha = 0.45) +
    geom_segment(data = dplyr::filter(cn_seg_abs, hl_abs),
                 aes(x = x_start, xend = x_end, y = cn, yend = cn, color = label),
                 linewidth = 1.8, alpha = 1) +
    scale_color_manual(values = cn_pal, name = "Timepoint") +
    scale_x_continuous(breaks = chrom_sizes$mid,
                       labels = str_remove(chrom_sizes$chromosome, "chr"),
                       expand = c(0.01, 0)) +
    scale_y_continuous(breaks = scales::breaks_width(1)) +
    coord_cartesian(ylim = c(0, cn_abs_cap)) +
    labs(
      x = "Chromosome",
      y = "Called copy number (cn)",
      title = paste0(clone, " called copy-number landscape across the ordered comparison series"),
      subtitle = paste0("Single panel, colored by timepoint | dashed = per-comparison ploidy baseline (",
                        paste(paste0(labels_ord, "=", label_ploidy[labels_ord]), collapse=", "),
                        ") | highlight = cn outside ploidy +/- 1 | y clipped at ", CN_CAP)
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x     = element_text(size = 7),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )

  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_cn_absolute.pdf")),
         p_cn_abs, width = 18, height = 6, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_manhattan_cn_absolute.png")),
         p_cn_abs, width = 18, height = 6, dpi = 150, limitsize = FALSE)
  message("absolute-CN Manhattan saved (y clipped at ", round(cn_abs_cap, 1), ")")
}

per_tp <- pmap(timeline, read_timepoint) %>% list_rbind()

if (nrow(per_tp) == 0) {
  message(clone, " temporal section: no called genes with cn in the timeline .cns files.")
  empty_gene_tp <- tibble(gene = character(), timepoint = integer(), label = character(),
                          chromosome = character(), start = integer(), end = integer(),
                          log2 = numeric(), cn = numeric(), ploidy = numeric(),
                          comparison = character(), reference = character(),
                          purity = numeric(), sample = character())
  empty_gene_pos <- tibble(gene = character(), chromosome = character(),
                           start = integer(), end = integer())
  empty_cn_wide <- tibble(gene = character())
  for (lbl in labels_ord) empty_cn_wide[[lbl]] <- numeric()
  drift_out <- tibble(gene = character(), cn_min = integer(), cn_max = integer(),
                      drift = integer(), n_tp = integer(), chromosome = character(),
                      start = integer(), end = integer(), is_oncogene = logical())
  for (lbl in labels_ord) {
    drift_out[[lbl]] <- integer()
    drift_out[[paste0("ploidy_", lbl)]] <- label_ploidy[lbl]
    drift_out[[paste0("purity_", lbl)]] <- label_purity[lbl]
  }
  write_csv(drift_out, file.path(TEMP_DIR, paste0(slug, "_gene_cn_drift.csv")))
  p_hm_genomic <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No called genes with cn", size = 4) +
    labs(title = paste0(clone, " copy-number drift across the ordered comparison series")) +
    theme_minimal(base_size = 10)
  p_hm_clust <- p_hm_genomic
  hm_w <- 8
  ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_genomic.pdf")), p_hm_genomic, width = hm_w, height = 4, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_genomic.png")), p_hm_genomic, width = hm_w, height = 4, dpi = 150, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_clustered.pdf")), p_hm_clust, width = hm_w, height = 4, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_clustered.png")), p_hm_clust, width = hm_w, height = 4, dpi = 150, limitsize = FALSE)
  p_traj <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No called genes with cn", size = 4) +
    labs(title = paste0(clone, " copy-number trajectories")) +
    theme_bw(base_size = 10)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_trajectories.pdf")), p_traj, width = 9, height = 6, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_trajectories.png")), p_traj, width = 9, height = 6, dpi = 150, limitsize = FALSE)
  p_allu <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No called genes with cn", size = 4) +
    labs(title = paste0(clone, " copy-number state flow across the ordered comparison series")) +
    theme_bw(base_size = 10)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_cn_state_alluvial.pdf")), p_allu, width = 9, height = 6, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_cn_state_alluvial.png")), p_allu, width = 9, height = 6, dpi = 150, limitsize = FALSE)
  panels <- Filter(Negate(is.null), list(p_cn_abs, p_hm_clust, p_hm_genomic, p_traj, p_allu))
  combined <- patchwork::wrap_plots(panels, ncol = 1)
  ggsave(file.path(OUT_DIR, paste0("Genomic_drift_", clone, ".pdf")),
         combined, width = 16, height = 4.5 * length(panels), limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0("Genomic_drift_", clone, ".png")),
         combined, width = 16, height = 4.5 * length(panels), dpi = 150, limitsize = FALSE)
  return(invisible(list(clone = clone, gene_tp = empty_gene_tp, gene_pos = empty_gene_pos,
                        cn_wide = empty_cn_wide, drift = drift_out,
                        labels_ord = labels_ord, ploidy = ploidy_of[[timeline$dir[1]]],
                        label_ploidy = label_ploidy)))
} else {

# Per gene per timepoint keep the segment CN furthest from neutral (per-comparison ploidy)
gene_tp <- per_tp %>%
  group_by(gene, timepoint, label) %>%
  slice_max(abs(cn - ploidy), n = 1, with_ties = FALSE) %>%
  ungroup()

# Representative genomic position per gene (chromosome with the largest span)
gene_pos <- per_tp %>%
  group_by(gene, chromosome) %>%
  summarise(start = min(start), end = max(end),
            span = sum(end - start), .groups = "drop_last") %>%
  slice_max(span, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, chromosome, start, end)

# Wide CN matrix (gene x timepoint); genes absent in a timepoint stay NA
cn_wide <- gene_tp %>%
  select(gene, label, cn) %>%
  tidyr::complete(gene, label = labels_ord) %>%
  pivot_wider(names_from = label, values_from = cn)

# Drift = span of CN across the timeline
drift <- gene_tp %>%
  group_by(gene) %>%
  summarise(cn_min = min(cn), cn_max = max(cn),
            drift = max(cn) - min(cn),
            n_tp  = n_distinct(timepoint), .groups = "drop") %>%
  arrange(desc(drift), desc(cn_max))

# Canonical oncogenes / lung-adeno drivers (defined globally above)
top_drift <- drift %>% slice_head(n = TOP_N) %>% pull(gene)
onco_present <- intersect(ONCOGENES, gene_tp$gene)
sel_genes <- union(top_drift, onco_present)

# Persist the drift table + wide CN matrix
drift_out <- drift %>% left_join(gene_pos, by = "gene") %>%
  left_join(cn_wide, by = "gene") %>%
  mutate(is_oncogene = gene %in% ONCOGENES) %>%
  arrange(desc(drift))
for (lbl in labels_ord) {
  drift_out[[paste0("ploidy_", lbl)]] <- label_ploidy[lbl]
  drift_out[[paste0("purity_", lbl)]] <- label_purity[lbl]
}
write_csv(drift_out, file.path(TEMP_DIR, paste0(slug, "_gene_cn_drift.csv")))
message(clone, " timeline genes: ", nrow(drift),
        " | selected for plots: ", length(sel_genes),
        " (top-drift ", length(top_drift), " + oncogenes ", length(onco_present), ")")

# -- 9a. Temporal heatmaps (gene x timepoint, fill = copy number) --------------
# Long data for the selected genes, CN clamped to [0,heatmap_cap] for the color scale.
hm_df <- gene_tp %>%
  filter(gene %in% sel_genes) %>%
  right_join(expand_grid(gene = sel_genes, label = labels_ord), by = c("gene","label")) %>%
  mutate(cn_cap = pmin(pmax(cn, 0), HEATMAP_CAP),
         label  = factor(label, levels = labels_ord),
         is_oncogene = gene %in% ONCOGENES)

# Gene ordering 1: by genomic position (chromosome, start)
gene_order_genomic <- gene_pos %>%
  filter(gene %in% sel_genes) %>%
  mutate(chr = factor(chromosome, levels = chrom_order)) %>%
  arrange(chr, start) %>% pull(gene)

# Gene ordering 2: hierarchical clustering of the CN profiles (complete rows only)
mat_df <- cn_wide %>% filter(gene %in% sel_genes)
mat <- as.matrix(mat_df[, labels_ord]); rownames(mat) <- mat_df$gene
gene_order_clust <- if (nrow(mat) > 2) {
  ok <- stats::complete.cases(mat)
  if (sum(ok) > 2) rownames(mat)[hclust(dist(mat[ok, , drop = FALSE]))$order]
  else gene_order_genomic
} else rownames(mat)

cn_heatmap <- function(gene_levels, subtitle) {
  ggplot(hm_df %>% mutate(gene = factor(gene, levels = gene_levels)),
         aes(x = gene, y = label, fill = cn_cap)) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = base_ploidy, limits = c(0, HEATMAP_CAP),
                         breaks = 0:HEATMAP_CAP, name = "Copy\nnumber") +
    scale_y_discrete(limits = rev(labels_ord)) +
    labs(x = NULL, y = NULL,
         title = paste0(clone, " copy-number drift across the ordered comparison series"),
         subtitle = subtitle) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6,
                                     face = ifelse(levels(factor(hm_df$gene, levels=gene_levels)) %in% ONCOGENES, "bold", "plain")),
          panel.grid = element_blank(),
          legend.position = "right")
}

p_hm_genomic <- cn_heatmap(gene_order_genomic,
  "Genes ordered by genomic coordinate | oncogenes in bold | top-drift + canonical drivers")
p_hm_clust <- cn_heatmap(gene_order_clust,
  "Genes ordered by hierarchical clustering of CN profiles | oncogenes in bold")

hm_w <- max(8, length(sel_genes) * 0.16)
ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_genomic.pdf")),  p_hm_genomic, width = hm_w, height = 4, limitsize = FALSE)
ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_genomic.png")),  p_hm_genomic, width = hm_w, height = 4, dpi = 150, limitsize = FALSE)
ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_clustered.pdf")), p_hm_clust,   width = hm_w, height = 4, limitsize = FALSE)
ggsave(file.path(TEMP_DIR, paste0(slug, "_temporal_heatmap_clustered.png")), p_hm_clust,   width = hm_w, height = 4, dpi = 150, limitsize = FALSE)

# -- 9b. Trajectory lines for the most dynamic genes + oncogenes ---------------
# A few tiny segments in repetitive regions (e.g. the chr7q22 MUC gene cluster)
# get spuriously high absolute CN (cn > 20-60) and otherwise flatten the axis.
# Clamp the plotted value so real cn 3-13 dynamics stay legible; clamped points
# sit at the CN_TRAJ_CAP ceiling.
# Needs >= 2 timepoints to draw a trajectory; single-timepoint clones skip this.
CN_TRAJ_CAP <- CN_CAP
p_traj <- NULL
if (length(labels_ord) >= 2) {
traj_genes <- union(head(top_drift, 15), onco_present)
traj_df <- hm_df %>%
  filter(gene %in% traj_genes) %>%
  mutate(tp = as.integer(factor(label, levels = labels_ord)),
         cn_plot = pmin(cn, CN_TRAJ_CAP))
traj_end <- traj_df %>% filter(tp == max(tp))
baseline_traj <- tibble(tp = seq_along(labels_ord), ploidy = label_ploidy[labels_ord])

p_traj <- ggplot(traj_df, aes(x = tp, y = cn_plot, group = gene)) +
  # background: all selected trajectories in grey
  geom_line(color = "grey75", linewidth = 0.4, alpha = 0.7) +
  geom_line(data = baseline_traj, aes(x = tp, y = ploidy, group = 1),
            linetype = "dashed", color = "black", linewidth = 0.3) +
  # foreground: oncogenes highlighted
  geom_line(data = filter(traj_df, is_oncogene),
            aes(color = gene), linewidth = 1) +
  geom_point(data = filter(traj_df, is_oncogene),
             aes(color = gene), size = 1.6) +
  {if (requireNamespace("ggrepel", quietly = TRUE))
     ggrepel::geom_text_repel(data = filter(traj_end, is_oncogene),
       aes(label = gene, color = gene), hjust = 0, size = 2.8,
       direction = "y", nudge_x = 0.15, segment.size = 0.2,
       max.overlaps = Inf, show.legend = FALSE)
   else
     geom_text(data = filter(traj_end, is_oncogene),
       aes(label = gene, color = gene), hjust = -0.15, size = 2.8, show.legend = FALSE)} +
  scale_x_continuous(breaks = seq_along(labels_ord), labels = labels_ord,
                     expand = expansion(mult = c(0.03, 0.30))) +
  coord_cartesian(ylim = c(0, CN_TRAJ_CAP)) +
  labs(x = NULL, y = "Copy number",
       title = paste0(clone, " copy-number trajectories"),
       subtitle = paste0("Grey = top-drift genes | colored = oncogenes | dashed = per-comparison ploidy baseline | y clamped at cn ",
                         CN_TRAJ_CAP, " (repetitive-region outliers sit at ceiling)"),
       color = "Oncogene") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(file.path(TEMP_DIR, paste0(slug, "_trajectories.pdf")), p_traj, width = 9, height = 6, limitsize = FALSE)
ggsave(file.path(TEMP_DIR, paste0(slug, "_trajectories.png")), p_traj, width = 9, height = 6, dpi = 150, limitsize = FALSE)
} else {
  message("Trajectories skipped for ", clone, " (needs >= 2 timepoints).")
}

# -- 9c. Alluvial of CN-state flow ---------------------------------------------
# Bins CN into states relative to each timepoint's own ploidy and shows how
# genes flow between states across timepoints; missing calls are "unknown".
p_allu <- NULL
if (length(labels_ord) >= 2) {
  state_lv <- c("loss", "neutral", "gain", "amplification", "unknown")
  cn_state <- function(x) cut(x, breaks = c(-Inf, -0.5, 0.5, 1.5, Inf),
                              labels = c("loss", "neutral", "gain", "amplification"))
  allu_df <- hm_df %>%
    mutate(ploidy_lbl = label_ploidy[as.character(label)],
           state = ifelse(is.na(cn), "unknown",
                          as.character(cn_state(cn - ploidy_lbl))),
           label = factor(label, levels = labels_ord))
  p_allu <- ggplot(allu_df,
                   aes(x = label, stratum = state, alluvium = gene,
                       fill = state, label = state)) +
    ggalluvial::geom_flow(alpha = 0.6) +
    ggalluvial::geom_stratum(alpha = 0.9) +
    scale_fill_manual(values = setNames(c("#2166AC","grey85","#F4A582","#B2182B","grey50"), state_lv)) +
    labs(x = NULL, y = "Number of genes", fill = "CN state",
          title = paste0(clone, " copy-number state flow across the ordered comparison series")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(TEMP_DIR, paste0(slug, "_cn_state_alluvial.pdf")), p_allu, width = 9, height = 6, limitsize = FALSE)
  ggsave(file.path(TEMP_DIR, paste0(slug, "_cn_state_alluvial.png")), p_allu, width = 9, height = 6, dpi = 150, limitsize = FALSE)
  message("Alluvial CN-state plot written.")
} else {
  message("Alluvial skipped for ", clone, " (needs >= 2 timepoints).")
}

# -- 9d. Combined multi-panel "genomic drift" figure ---------------------------
# Stack the per-clone temporal panels (absolute-CN Manhattan, both drift
# heatmaps, trajectories, CN-state alluvial) into a single overview figure.
panels <- Filter(Negate(is.null), list(p_cn_abs, p_hm_clust, p_hm_genomic, p_traj, p_allu))
if (length(panels) > 0) {
  combined <- patchwork::wrap_plots(panels, ncol = 1)
  ggsave(file.path(OUT_DIR, paste0("Genomic_drift_", clone, ".pdf")),
         combined, width = 16, height = 4.5 * length(panels), limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0("Genomic_drift_", clone, ".png")),
         combined, width = 16, height = 4.5 * length(panels), dpi = 150, limitsize = FALSE)
  message(clone, " combined genomic-drift figure written to ", OUT_DIR)
}

message(clone, " temporal outputs written to ", TEMP_DIR)

invisible(list(clone = clone, gene_tp = gene_tp, gene_pos = gene_pos,
               cn_wide = cn_wide, drift = drift, labels_ord = labels_ord,
               ploidy = ploidy_of[[timeline$dir[1]]], label_ploidy = label_ploidy))

}  # end per-clone temporal body (per_tp non-empty)
}  # end run_clone_temporal()

# Run each configured clone timeline; collect results in a list keyed by clone.
temporal_results <- list()
for (spec in TEMPORAL_CLONES) {
  temporal_results[[spec$clone]] <-
    run_clone_temporal(spec$clone, spec$timeline, spec$out_subdir)
}
if (length(TEMPORAL_CLONES) == 0)
  message("Section 9: no TEMPORAL_CLONES configured - skipping temporal analysis.")


# ==============================================================================
# 10. Conserved copy-number alterations shared by a pair of clones
# ==============================================================================
# Both clones are compared against the same human reference, so their vs-reference
# integer CN is directly comparable. Gains/losses are called relative to each
# clone's OWN modal ploidy (near-triploid) so a whole-genome baseline shift is
# not mistaken for gene-level change. A gene is "conserved" when both clones move
# in the same direction (both gained or both lost). Driven by CONSERVED_COMPARISONS.
run_conserved_pair <- function(res_a, res_b, name_a, name_b, out_subdir) {
  cons_dir <- file.path(OUT_DIR, out_subdir)
  dir.create(cons_dir, recursive = TRUE, showWarnings = FALSE)
  slug <- tolower(paste0(name_a, "_", name_b))

  empty_cons <- tibble(gene = character(),
                       !!paste0("cn_", name_a) := integer(),
                       !!paste0("cn_", name_b) := integer(),
                       !!paste0("state_", name_a) := character(),
                       !!paste0("state_", name_b) := character(),
                       conserved_state = character(), chromosome = character(),
                       start = integer(), end = integer(), is_oncogene = logical(),
                       chr = factor())
  write_empty <- function() {
    write_csv(empty_cons, file.path(cons_dir, paste0(slug, "_gene_cn_all.csv")))
    write_csv(empty_cons, file.path(cons_dir, paste0(slug, "_conserved_cnv.csv")))
    p_sc <- ggplot() +
      annotate("text", x = CN_CAP/2, y = CN_CAP/2, label = "No called genes with cn", size = 4) +
      coord_cartesian(xlim = c(0, CN_CAP), ylim = c(0, CN_CAP)) +
      labs(x = paste0(name_a, " called CN (vs reference)"),
           y = paste0(name_b, " called CN (vs reference)"),
           title = paste0("Conserved called CN: ", name_a, " vs ", name_b)) +
      theme_bw(base_size = 10)
    ggsave(file.path(cons_dir, paste0(slug, "_cn_scatter.pdf")), p_sc, width = 7, height = 6, limitsize = FALSE)
    ggsave(file.path(cons_dir, paste0(slug, "_cn_scatter.png")), p_sc, width = 7, height = 6, dpi = 150, limitsize = FALSE)
    save_scatter_html(empty_cons %>% mutate(x = .data[[paste0("cn_", name_a)]],
                                            y = .data[[paste0("cn_", name_b)]],
                                            state = conserved_state),
                      paste0(name_a, " called CN (vs reference)"),
                      paste0(name_b, " called CN (vs reference)"),
                      paste0("Conserved called CN: ", name_a, " vs ", name_b),
                      file.path(cons_dir, paste0(slug, "_cn_scatter")),
                      c("conserved gain" = "#B2182B"), cap = CN_CAP)
    p_hm <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No called genes with cn", size = 4) +
      labs(title = paste0("Conserved CNV genes shared by ", name_a, " and ", name_b)) +
      theme_minimal(base_size = 9)
    ggsave(file.path(cons_dir, paste0(slug, "_conserved_heatmap.pdf")), p_hm, width = 8, height = 3, limitsize = FALSE)
    ggsave(file.path(cons_dir, paste0(slug, "_conserved_heatmap.png")), p_hm, width = 8, height = 3, dpi = 150, limitsize = FALSE)
    message("Conserved ", name_a, "/", name_b, " outputs written to ", cons_dir)
  }

  if (is.null(res_a) || is.null(res_b) || is.null(res_a$gene_tp) || is.null(res_b$gene_tp) ||
      nrow(res_a$gene_tp) == 0 || nrow(res_b$gene_tp) == 0) {
    write_empty()
    return(invisible(NULL))
  }

  # Per-gene called CN at each clone's reference (first) timepoint.
  ref_cn <- function(res) {
    res$gene_tp %>%
      filter(label == res$labels_ord[1]) %>%
      group_by(gene) %>% slice_max(abs(cn - ploidy), n = 1, with_ties = FALSE) %>%
      ungroup() %>% select(gene, cn)
  }
  a <- ref_cn(res_a) %>% rename(cn_a = cn)
  b <- ref_cn(res_b) %>% rename(cn_b = cn)
  base_a <- res_a$ploidy
  base_b <- res_b$ploidy

  gene_pos_all <- bind_rows(res_a$gene_pos, res_b$gene_pos) %>%
    distinct(gene, .keep_all = TRUE)

  call_state <- function(cn, base) case_when(cn >= base + 1 ~ "gain",
                                             cn <= base - 1 ~ "loss",
                                             TRUE            ~ "neutral")
  conserved <- inner_join(a, b, by = "gene") %>%
    mutate(state_a = call_state(cn_a, base_a),
           state_b = call_state(cn_b, base_b),
           conserved_state = case_when(
             state_a == "gain"    & state_b == "gain"    ~ "conserved gain",
             state_a == "loss"    & state_b == "loss"    ~ "conserved loss",
             state_a == "neutral" & state_b == "neutral" ~ "neutral in both",
             TRUE                                        ~ "divergent")) %>%
    left_join(gene_pos_all, by = "gene") %>%
    mutate(is_oncogene = gene %in% ONCOGENES,
           chr = factor(chromosome, levels = chrom_order)) %>%
    arrange(chr, start)

  if (nrow(conserved) == 0) {
    write_empty()
    return(invisible(NULL))
  }

  # Rename the generic cn_/state_ columns to the clone names for output.
  out_names <- c(cn_a = paste0("cn_", name_a),   cn_b = paste0("cn_", name_b),
                 state_a = paste0("state_", name_a), state_b = paste0("state_", name_b))
  conserved_out <- conserved %>% rename(!!!setNames(names(out_names), out_names))
  write_csv(conserved_out, file.path(cons_dir, paste0(slug, "_gene_cn_all.csv")))
  cons_alt <- conserved %>%
    filter(conserved_state %in% c("conserved gain", "conserved loss")) %>%
    arrange(conserved_state, chr, start)
  cons_alt %>% rename(!!!setNames(names(out_names), out_names)) %>%
    write_csv(file.path(cons_dir, paste0(slug, "_conserved_cnv.csv")))
  message("Conserved CNV (", name_a, "/", name_b, "): ",
          sum(conserved$conserved_state == "conserved gain"), " gains, ",
          sum(conserved$conserved_state == "conserved loss"),
          " losses (configured ploidy ", name_a, "=", base_a, ", ", name_b, "=", base_b, ")")

  # -- 10a. Per-gene CN scatter: name_a vs name_b (colored by conserved state) --
  cap <- CN_CAP
  sc_df <- conserved %>% mutate(x = pmin(cn_a, cap), y = pmin(cn_b, cap))
  state_cols <- c("conserved gain" = "#B2182B", "conserved loss" = "#2166AC",
                  "divergent" = "grey60", "neutral in both" = "grey85")
  p_sc <- ggplot(sc_df, aes(x, y, color = conserved_state)) +
    geom_count(alpha = 0.6) +
    scale_size_area(max_size = 9, name = "genes") +
    geom_vline(xintercept = base_a, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = base_b, linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = state_cols, name = "State") +
    {if (requireNamespace("ggrepel", quietly = TRUE))
       ggrepel::geom_text_repel(
         data = filter(sc_df, is_oncogene &
                       conserved_state %in% c("conserved gain","conserved loss")),
         aes(label = gene), size = 2.6, color = "black", fontface = "bold",
         max.overlaps = Inf, show.legend = FALSE, inherit.aes = TRUE,
         min.segment.length = 0, segment.size = 0.3, segment.color = "grey40",
         box.padding = 0.6, point.padding = 0.3, force = 3, seed = 1,
         bg.color = "white", bg.r = 0.18)} +
    coord_cartesian(xlim = c(0, cap), ylim = c(0, cap)) +
    labs(x = paste0(name_a, " called CN (vs reference)"),
         y = paste0(name_b, " called CN (vs reference)"),
         title = paste0("Conserved called CN: ", name_a, " vs ", name_b),
         subtitle = paste0("Dashed = configured ploidy (first reference comparison) | axes clamped at ", cap,
                           " | conserved-altered oncogenes labelled")) +
    theme_bw(base_size = 10)
  ggsave(file.path(cons_dir, paste0(slug, "_cn_scatter.pdf")), p_sc, width = 7, height = 6, limitsize = FALSE)
  ggsave(file.path(cons_dir, paste0(slug, "_cn_scatter.png")), p_sc, width = 7, height = 6, dpi = 150, limitsize = FALSE)
  save_scatter_html(conserved %>% mutate(x = cn_a, y = cn_b, state = conserved_state),
                    paste0(name_a, " called CN (vs reference)"),
                    paste0(name_b, " called CN (vs reference)"),
                    paste0("Conserved called CN: ", name_a, " vs ", name_b),
                    file.path(cons_dir, paste0(slug, "_cn_scatter")), state_cols, cap = cap)

  # -- 10b. Heatmap of the conserved-altered genes (both clones x genes) --------
  # The full conserved list can run to thousands of genes (both lines are near-
  # triploid), which is neither readable nor renderable as a single heatmap. Show
  # a curated subset: all conserved-altered oncogenes plus the strongest deviations
  # by combined distance from each clone's modal ploidy. Full lists stay in the CSVs.
  if (nrow(cons_alt) > 0) {
    cons_rank <- cons_alt %>%
      mutate(dev = abs(cn_a - base_a) + abs(cn_b - base_b))
    hm_genes <- cons_rank %>%
      filter(is_oncogene) %>%
      bind_rows(cons_rank %>% filter(!is_oncogene) %>% slice_max(dev, n = HM_MAX, with_ties = FALSE)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      arrange(chr, start)
    hm_cons <- hm_genes %>%
      transmute(gene, !!name_a := cn_a, !!name_b := cn_b) %>%
      tidyr::pivot_longer(c(!!name_a, !!name_b), names_to = "clone", values_to = "cn") %>%
      mutate(cn_cap = pmin(pmax(cn, 0), HEATMAP_CAP),
             gene = factor(gene, levels = unique(hm_genes$gene)))
    p_hm <- ggplot(hm_cons, aes(x = gene, y = clone, fill = cn_cap)) +
      geom_tile(color = "grey90", linewidth = 0.2) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = round((base_a + base_b) / 2), limits = c(0, HEATMAP_CAP),
                           breaks = 0:HEATMAP_CAP, name = "Copy\nnumber") +
      labs(x = NULL, y = NULL,
           title = paste0("Conserved CNV genes shared by ", name_a, " and ", name_b),
           subtitle = paste0("Called CN (vs reference), capped at ", HEATMAP_CAP, " | oncogenes + top ",
                             HM_MAX, " by deviation | full lists in CSV")) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6,
              face = ifelse(levels(hm_cons$gene) %in% ONCOGENES, "bold", "plain")),
            panel.grid = element_blank())
    hm_w <- max(8, nlevels(hm_cons$gene) * 0.16)
    ggsave(file.path(cons_dir, paste0(slug, "_conserved_heatmap.pdf")), p_hm, width = hm_w, height = 3, limitsize = FALSE)
    ggsave(file.path(cons_dir, paste0(slug, "_conserved_heatmap.png")), p_hm, width = hm_w, height = 3, dpi = 150, limitsize = FALSE)
  } else {
    p_hm <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No conserved altered genes", size = 4) +
      labs(title = paste0("Conserved CNV genes shared by ", name_a, " and ", name_b)) +
      theme_minimal(base_size = 9)
    ggsave(file.path(cons_dir, paste0(slug, "_conserved_heatmap.pdf")), p_hm, width = 8, height = 3, limitsize = FALSE)
    ggsave(file.path(cons_dir, paste0(slug, "_conserved_heatmap.png")), p_hm, width = 8, height = 3, dpi = 150, limitsize = FALSE)
  }
  message("Conserved ", name_a, "/", name_b, " outputs written to ", cons_dir)
}

for (cmp in CONSERVED_COMPARISONS) {
  ra <- temporal_results[[cmp$clones[1]]]
  rb <- temporal_results[[cmp$clones[2]]]
  run_conserved_pair(ra, rb, cmp$clones[1], cmp$clones[2], cmp$out_subdir)
}


# ==============================================================================
# 11. Additional per-gene CN scatters (same style as the conserved-CNV plot)
# ==============================================================================
# Each entry in SCATTER_COMPARISONS declares two comparison directories that
# provide the x/y axes; per-gene integer CN is read from each and plotted with
# cn_state_scatter(). Add/remove entries in the USER CONFIGURATION block.
scatter_dir <- file.path(OUT_DIR, "cn_scatters")
dir.create(scatter_dir, showWarnings = FALSE)

# Collapse per-gene CN (furthest from neutral) for a two-column comparison and
# hand it to cn_state_scatter(). 'gene_tp' has columns gene,label,cn (plus
# chromosome,start used to annotate the interactive hover).
scatter_from_gene_tp <- function(gene_tp, xlab_lbl, ylab_lbl, xlab, ylab, title, out_prefix,
                                 bx = NULL, by = NULL) {
  gene_pos <- gene_tp %>% distinct(gene, .keep_all = TRUE) %>%
    transmute(gene, chr = chromosome, start)
  tab <- gene_tp %>%
    filter(label %in% c(xlab_lbl, ylab_lbl)) %>%
    group_by(gene, label) %>% slice_max(abs(cn - ploidy), n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(gene, label, cn) %>%
    pivot_wider(names_from = label, values_from = cn)
  for (lbl in c(xlab_lbl, ylab_lbl)) if (!lbl %in% names(tab)) tab[[lbl]] <- NA_real_
  tab <- tab %>%
    filter(!is.na(.data[[xlab_lbl]]), !is.na(.data[[ylab_lbl]])) %>%
    transmute(gene, x = .data[[xlab_lbl]], y = .data[[ylab_lbl]],
              is_oncogene = gene %in% ONCOGENES) %>%
    left_join(gene_pos, by = "gene")
  cn_state_scatter(tab, xlab, ylab, title, out_prefix, bx = bx, by = by)
  message("Scatter written: ", out_prefix, ".png (", nrow(tab), " genes)")
}

# Build a scatter table from two genemetrics TSVs (called, ploidy-corrected CN).
# Returns a tibble compatible with cn_state_scatter().
scatter_from_genometrics <- function(path_x, path_y, oncogenes, ploidy_x, ploidy_y) {
  if (!file.exists(path_x)) stop("Genemetrics file not found: ", path_x)
  if (!file.exists(path_y)) stop("Genemetrics file not found: ", path_y)
  read_gm <- function(p, ploidy) {
    read_typed_tsv(p, required = c("gene","chromosome","start","end","log2","weight","cn"),
                   numeric_cols = c("start","end","probes","log2","weight","cn")) %>%
      transmute(gene = gene,
                chr  = chromosome,
                start = as.integer(start),
                cn   = as.numeric(cn)) %>%
      filter(!is.na(cn)) %>%
      group_by(gene) %>% slice_max(abs(cn - ploidy), n = 1, with_ties = FALSE) %>% ungroup()
  }
  gx <- read_gm(path_x, ploidy_x)
  gy <- read_gm(path_y, ploidy_y)
  tab <- inner_join(gx, gy, by = "gene", suffix = c("_x", "_y")) %>%
    transmute(gene,
              chr   = chr_x,
              start = start_x,
              x     = cn_x,
              y     = cn_y,
              is_oncogene = gene %in% oncogenes)
  tab
}

# -- Run each configured scatter (two comparisons -> x/y axes) -----------------
for (sc in SCATTER_COMPARISONS) {
  tl <- sc$timeline
  if (!"timepoint" %in% names(tl)) tl$timepoint <- seq_len(nrow(tl))
  tp <- pmap(tl[c("timepoint", "label", "dir")], read_timepoint) %>% list_rbind()
  scatter_from_gene_tp(tp, tl$label[1], tl$label[2], sc$xlab, sc$ylab, sc$title,
                       file.path(scatter_dir, sc$out_prefix),
                       bx = ploidy_of[[tl$dir[1]]], by = ploidy_of[[tl$dir[2]]])
  if (!is.null(sc$gm)) {
    gm_xf <- gm_of[[sc$gm$x]]
    gm_yf <- gm_of[[sc$gm$y]]
    if (is.null(gm_xf) || is.na(gm_xf) || is.null(gm_yf) || is.na(gm_yf)) {
      stop("Absolute scatter: unknown gm comparison: ", sc$out_prefix)
    }
    tab_abs <- scatter_from_genometrics(gm_xf, gm_yf, ONCOGENES,
                                        ploidy_of[[sc$gm$x]], ploidy_of[[sc$gm$y]])
    if (identical(sc$gm$x, tl$dir[1]) && identical(sc$gm$y, tl$dir[2])) {
      abs_xlab  <- gsub("\\s*\\(vs[^)]*\\)", " (called CN)", sc$xlab)
      abs_ylab  <- gsub("\\s*\\(vs[^)]*\\)", " (called CN)", sc$ylab)
    } else {
      abs_xlab  <- paste0(sc$gm$x, " called CN")
      abs_ylab  <- paste0(sc$gm$y, " called CN")
    }
    abs_title <- paste0(sc$title, " - called CN")
    cn_state_scatter(tab_abs, abs_xlab, abs_ylab, abs_title,
                     file.path(scatter_dir, paste0(sc$out_prefix, "_abs")), cap = CN_CAP,
                     bx = ploidy_of[[sc$gm$x]], by = ploidy_of[[sc$gm$y]])
    message("Absolute scatter written: ", sc$out_prefix, "_abs.png (", nrow(tab_abs), " genes)")
  }
}
if (length(SCATTER_COMPARISONS) == 0)
  message("Section 11: no SCATTER_COMPARISONS configured - skipping scatters.")

sizes <- file.info(expected_outputs)$size
missing_outputs <- expected_outputs[is.na(sizes) | sizes == 0]
if (length(missing_outputs) > 0)
  stop("Expected outputs not generated or empty: ", paste(missing_outputs, collapse=", "))
message("All ", length(expected_outputs), " expected outputs present.")


