# ============================================================
# Fig. 4B — V-associated bulk–rhizosphere bacterial similarity
# ============================================================
#
# HOW TO USE
# 1) Put this .R file in:
#    C:/Users/lenovo/OneDrive/01 V-ARG submission20260527/01 Data analysis/fig 4b
# 2) Keep the two source files somewhere under the same folder:
#      - Soil geochemistry.xlsx
#      - *asv.full.xls   (e.g. ASV_20240225_195022_asv.full.xls)
# 3) Open this script in RStudio and Run All.
#
# OUTPUT
# All results are written to:
#    ./output/
#
# Main outputs:
#    Fig4B_source_data.csv
#    Fig4B_statistics.csv
#    Fig4B_model_summary.txt
#    Fig4B_similarity_vs_V.pdf
#    Fig4B_similarity_vs_V.svg
#    Fig4B_similarity_vs_V.png
#
# ANALYSIS
# - Pairing:
#     A1-x ↔ A2-x
#     B1-x ↔ B2-x
#     C1-x ↔ C2-x
#   where "1" = rhizosphere and "2" = bulk soil.
#
# - Bray-Curtis similarity = 1 - Bray-Curtis distance
# - Jaccard similarity     = 1 - binary Jaccard distance
#
# - Exposure:
#     log10(bulk-soil V)
#
# - Adjusted model:
#     standardized similarity ~ standardized log10(V)
#                               + region
#                               + standardized bulk-soil pH
#                               + standardized bulk-soil AP
#
# - P value:
#     Freedman-Lane permutation test (9,999 permutations),
#     permuting reduced-model residuals within region.
#
# - Figure:
#     raw pairwise observations + separate lm fitted lines / 95% CI.
#     The fitted lines are visual summaries; the annotated beta and P
#     come from the adjusted permutation models above.
#
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "data.table",
  "readxl",
  "dplyr",
  "tidyr",
  "vegan",
  "ggplot2",
  "svglite"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  message(
    "Installing missing R packages: ",
    paste(missing_pkgs, collapse = ", ")
  )

  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org",
    dependencies = NA
  )

  still_missing <- missing_pkgs[
    !vapply(missing_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(still_missing) > 0) {
    stop(
      "Package installation did not complete successfully: ",
      paste(still_missing, collapse = ", "),
      "\nInstall these packages in the RStudio Console and run the script again."
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
  library(svglite)
})

set.seed(20260915)


# ------------------------------------------------------------
# 1. Locate the folder containing this R script
# ------------------------------------------------------------

get_script_dir <- function() {

  # Preferred route when running inside RStudio
  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()
  ) {
    p <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )

    if (!is.null(p) && nzchar(p)) {
      return(dirname(normalizePath(p, winslash = "/", mustWork = TRUE)))
    }
  }

  # Fallback
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

BASE_DIR <- get_script_dir()

cat("\nWorking analysis folder:\n", BASE_DIR, "\n")


# ------------------------------------------------------------
# 2. Locate source files under BASE_DIR
# ------------------------------------------------------------

geo_hits <- list.files(
  path = BASE_DIR,
  pattern = "^Soil geochemistry\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(geo_hits) == 0) {
  # fallback in case the workbook name is slightly different
  geo_hits <- list.files(
    path = BASE_DIR,
    pattern = "geochemistry.*\\.xlsx$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

asv_hits <- list.files(
  path = BASE_DIR,
  pattern = "asv\\.full\\.xls$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(geo_hits) == 0) {
  stop(
    "Cannot find the soil geochemistry workbook under:\n",
    BASE_DIR
  )
}

if (length(asv_hits) == 0) {
  stop(
    "Cannot find the ASV table (*asv.full.xls) under:\n",
    BASE_DIR
  )
}

GEOCHEM_FILE <- geo_hits[1]
ASV_FILE     <- asv_hits[1]

cat("\nGeochemistry source:\n", GEOCHEM_FILE, "\n")
cat("\nASV source:\n", ASV_FILE, "\n")


# ------------------------------------------------------------
# 3. Output folder
# ------------------------------------------------------------

OUT_DIR <- file.path(BASE_DIR, "output")

if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

OUT_SOURCE <- file.path(OUT_DIR, "Fig4B_source_data.csv")
OUT_STATS  <- file.path(OUT_DIR, "Fig4B_statistics.csv")
OUT_TXT    <- file.path(OUT_DIR, "Fig4B_model_summary.txt")

OUT_PDF <- file.path(OUT_DIR, "Fig4B_similarity_vs_V.pdf")
OUT_SVG <- file.path(OUT_DIR, "Fig4B_similarity_vs_V.svg")
OUT_PNG <- file.path(OUT_DIR, "Fig4B_similarity_vs_V.png")


# ------------------------------------------------------------
# 4. Figure 4a-derived color family
# ------------------------------------------------------------
#
# Colors sampled from the supplied Fig. 4a:
#   Bulk soil    = #DBB999
#   Rhizosphere  = #F6EDC1
#   Endophyte    = #D9B6C9
#   Plant surface= #CAEDB4
#
# For Fig. 4B we use darker readable derivatives for points/lines,
# while retaining the original pastel colors for confidence ribbons.

COL_BRAY_LINE <- "#B27697"   # muted mauve
COL_BRAY_FILL <- "#D9B6C9"   # Fig.4a endophyte pastel

COL_JAC_LINE  <- "#78A85F"   # muted green
COL_JAC_FILL  <- "#CAEDB4"   # Fig.4a plant-surface pastel


# ------------------------------------------------------------
# 5. Utility functions
# ------------------------------------------------------------

norm_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("_", "-", x, fixed = TRUE)
  toupper(x)
}

region_from_id <- function(x) {
  dplyr::recode(
    substr(x, 1, 1),
    "A" = "Chengde",
    "B" = "Huaihua",
    "C" = "Panzhihua",
    .default = NA_character_
  )
}

read_asv_file <- function(path) {

  # Most Majorbio .xls ASV tables are tab-delimited text despite the extension.
  ans <- tryCatch(
    data.table::fread(
      file = path,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      data.table = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(ans) && ncol(ans) > 2) {
    return(ans)
  }

  # Fallback if it is a genuine Excel workbook
  ans <- tryCatch(
    as.data.frame(
      readxl::read_excel(
        path,
        .name_repair = "minimal"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(ans)) {
    stop("Unable to read ASV file:\n", path)
  }

  ans
}

pair_distance <- function(mat, id1, id2, method, binary = FALSE) {

  if (!all(c(id1, id2) %in% rownames(mat))) {
    return(NA_real_)
  }

  z <- mat[c(id1, id2), , drop = FALSE]

  if (any(rowSums(z) == 0)) {
    return(NA_real_)
  }

  as.numeric(
    vegan::vegdist(
      z,
      method = method,
      binary = binary
    )
  )
}


# ------------------------------------------------------------
# 6. Read and clean soil geochemistry
# ------------------------------------------------------------

geo_raw <- as.data.frame(
  readxl::read_excel(
    GEOCHEM_FILE,
    sheet = 1,
    .name_repair = "minimal"
  )
)

# exact columns confirmed from the supplied workbook:
# ID, PH, AP, ..., V, ...
names(geo_raw) <- trimws(names(geo_raw))

required_geo <- c("ID", "PH", "AP", "V")

if (!all(required_geo %in% names(geo_raw))) {
  stop(
    "Required columns not found in Soil geochemistry.xlsx.\n",
    "Expected: ID, PH, AP, V\n",
    "Found: ", paste(names(geo_raw), collapse = ", ")
  )
}

geo <- geo_raw %>%
  transmute(
    sample_id = norm_id(ID),
    pH = suppressWarnings(as.numeric(PH)),
    AP = suppressWarnings(as.numeric(AP)),
    V  = suppressWarnings(as.numeric(V))
  ) %>%
  filter(
    grepl("^[ABC][12]-[0-9]+$", sample_id)
  ) %>%
  distinct(sample_id, .keep_all = TRUE)

cat("\nSoil geochemistry samples retained:", nrow(geo), "\n")


# ------------------------------------------------------------
# 7. Read ASV table and extract soil samples
# ------------------------------------------------------------

asv_raw <- read_asv_file(ASV_FILE)

asv_names <- names(asv_raw)
asv_ids   <- norm_id(asv_names)

soil_col_idx <- which(
  grepl("^[ABC][12]-[0-9]+$", asv_ids)
)

if (length(soil_col_idx) == 0) {
  stop(
    "No soil sample columns (A1/A2/B1/B2/C1/C2) were detected ",
    "in the ASV table."
  )
}

asv_mat <- as.matrix(
  asv_raw[, soil_col_idx, drop = FALSE]
)

suppressWarnings(
  storage.mode(asv_mat) <- "numeric"
)

asv_mat[is.na(asv_mat)] <- 0

# transpose to samples × ASVs
asv_mat <- t(asv_mat)

rownames(asv_mat) <- asv_ids[soil_col_idx]

# remove all-zero ASVs
asv_mat <- asv_mat[
  ,
  colSums(asv_mat) > 0,
  drop = FALSE
]

# remove all-zero samples, if any
asv_mat <- asv_mat[
  rowSums(asv_mat) > 0,
  ,
  drop = FALSE
]

cat("Soil ASV samples retained:", nrow(asv_mat), "\n")
cat("ASVs retained:", ncol(asv_mat), "\n")


# ------------------------------------------------------------
# 8. Relative-abundance matrix for Bray-Curtis
# ------------------------------------------------------------

asv_rel <- sweep(
  asv_mat,
  1,
  rowSums(asv_mat),
  "/"
)


# ------------------------------------------------------------
# 9. Construct matched bulk–rhizosphere pairs
# ------------------------------------------------------------
#
# Pairing convention:
# A1-x ↔ A2-x
# B1-x ↔ B2-x
# C1-x ↔ C2-x
#
# 1 = rhizosphere
# 2 = bulk soil

available_ids <- intersect(
  rownames(asv_mat),
  geo$sample_id
)

pair_map <- data.frame(
  sample_id = available_ids,
  stringsAsFactors = FALSE
) %>%
  mutate(
    region_code = substr(sample_id, 1, 1),
    soil_type   = substr(sample_id, 2, 2),
    pair_no = as.integer(
      sub("^[ABC][12]-", "", sample_id)
    )
  ) %>%
  select(
    region_code,
    pair_no,
    soil_type,
    sample_id
  ) %>%
  pivot_wider(
    names_from = soil_type,
    values_from = sample_id,
    names_prefix = "type_"
  ) %>%
  filter(
    !is.na(type_1),
    !is.na(type_2)
  ) %>%
  transmute(
    pair_id = paste0(region_code, "-", pair_no),
    region = recode(
      region_code,
      "A" = "Chengde",
      "B" = "Huaihua",
      "C" = "Panzhihua"
    ),
    rhizosphere_sample = type_1,
    bulk_sample = type_2
  ) %>%
  arrange(region, pair_id)

cat("\nMatched bulk–rhizosphere pairs detected:", nrow(pair_map), "\n")


# ------------------------------------------------------------
# 10. Calculate Bray-Curtis and Jaccard similarities
# ------------------------------------------------------------

source_data <- pair_map %>%

  left_join(
    geo %>%
      rename(
        bulk_sample = sample_id,
        bulk_pH = pH,
        bulk_AP = AP,
        bulk_V = V
      ),
    by = "bulk_sample"
  ) %>%

  left_join(
    geo %>%
      rename(
        rhizosphere_sample = sample_id,
        rhizo_pH = pH,
        rhizo_AP = AP,
        rhizo_V = V
      ),
    by = "rhizosphere_sample"
  ) %>%

  rowwise() %>%

  mutate(

    Bray_distance = pair_distance(
      asv_rel,
      rhizosphere_sample,
      bulk_sample,
      method = "bray",
      binary = FALSE
    ),

    Jaccard_distance = pair_distance(
      asv_mat,
      rhizosphere_sample,
      bulk_sample,
      method = "jaccard",
      binary = TRUE
    ),

    Bray_similarity =
      1 - Bray_distance,

    Jaccard_similarity =
      1 - Jaccard_distance,

    log10_bulk_V =
      ifelse(
        !is.na(bulk_V) & bulk_V > 0,
        log10(bulk_V),
        NA_real_
      )
  ) %>%

  ungroup()

write.csv(
  source_data,
  OUT_SOURCE,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "Complete matched pairs with both similarities:",
  sum(
    complete.cases(
      source_data[, c(
        "Bray_similarity",
        "Jaccard_similarity"
      )]
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 11. Adjusted Freedman-Lane permutation test
# ------------------------------------------------------------

freedman_lane <- function(dat, response, nperm = 9999) {

  d <- dat %>%
    transmute(
      Y = .data[[response]],
      logV = log10_bulk_V,
      region = factor(region),
      pH = bulk_pH,
      AP = bulk_AP
    ) %>%
    filter(
      complete.cases(.)
    )

  if (nrow(d) < 10) {
    stop("Too few complete pairs for adjusted analysis of ", response)
  }

  # Fully standardized continuous variables
  d <- d %>%
    mutate(
      Y_z    = as.numeric(scale(Y)),
      logV_z = as.numeric(scale(logV)),
      pH_z   = as.numeric(scale(pH)),
      AP_z   = as.numeric(scale(AP))
    )

  full_fit <- lm(
    Y_z ~ logV_z + region + pH_z + AP_z,
    data = d
  )

  reduced_fit <- lm(
    Y_z ~ region + pH_z + AP_z,
    data = d
  )

  beta_obs <- unname(
    coef(full_fit)["logV_z"]
  )

  fitted0 <- fitted(reduced_fit)
  resid0  <- residuals(reduced_fit)

  beta_perm <- numeric(nperm)

  # Region-restricted permutation of reduced-model residuals
  region_index <- split(
    seq_len(nrow(d)),
    d$region
  )

  for (i in seq_len(nperm)) {

    perm_resid <- resid0

    for (idx in region_index) {
      perm_resid[idx] <- sample(
        resid0[idx],
        size = length(idx),
        replace = FALSE
      )
    }

    y_perm <- fitted0 + perm_resid

    fit_perm <- lm(
      y_perm ~ logV_z + region + pH_z + AP_z,
      data = d
    )

    beta_perm[i] <- unname(
      coef(fit_perm)["logV_z"]
    )
  }

  p_perm <- (
    sum(abs(beta_perm) >= abs(beta_obs)) + 1
  ) / (nperm + 1)

  sm <- summary(full_fit)

  list(
    N = nrow(d),
    beta = beta_obs,
    permutation_P = p_perm,
    parametric_P = coef(sm)["logV_z", "Pr(>|t|)"],
    adjusted_R2 = sm$adj.r.squared,
    model = full_fit,
    beta_perm = beta_perm
  )
}


res_bray <- freedman_lane(
  source_data,
  "Bray_similarity",
  nperm = 9999
)

res_jac <- freedman_lane(
  source_data,
  "Jaccard_similarity",
  nperm = 9999
)


stats_out <- data.frame(
  Metric = c(
    "Bray-Curtis similarity",
    "Jaccard similarity"
  ),
  N = c(
    res_bray$N,
    res_jac$N
  ),
  Standardized_beta_log10V = c(
    res_bray$beta,
    res_jac$beta
  ),
  Permutation_P = c(
    res_bray$permutation_P,
    res_jac$permutation_P
  ),
  Parametric_P = c(
    res_bray$parametric_P,
    res_jac$parametric_P
  ),
  Adjusted_R2 = c(
    res_bray$adjusted_R2,
    res_jac$adjusted_R2
  )
)

write.csv(
  stats_out,
  OUT_STATS,
  row.names = FALSE
)

cat("\nAdjusted analysis results:\n")
print(stats_out)


# ------------------------------------------------------------
# 12. Save model summaries
# ------------------------------------------------------------

sink(OUT_TXT)

cat("FIGURE 4B — MODEL SUMMARY\n")
cat("========================================\n\n")

cat("Analysis model:\n")
cat(
  "standardized similarity ~ standardized log10(bulk-soil V) + ",
  "region + standardized bulk-soil pH + standardized bulk-soil AP\n\n"
)

cat("Permutation:\n")
cat(
  "Freedman-Lane residual permutation; 9,999 permutations; ",
  "residuals permuted within region.\n\n"
)

cat("----------------------------------------\n")
cat("Bray-Curtis similarity\n")
cat("----------------------------------------\n")
print(summary(res_bray$model))
cat("\nPermutation P =", res_bray$permutation_P, "\n\n")

cat("----------------------------------------\n")
cat("Jaccard similarity\n")
cat("----------------------------------------\n")
print(summary(res_jac$model))
cat("\nPermutation P =", res_jac$permutation_P, "\n")

sink()


# ------------------------------------------------------------
# 13. Prepare plotting data
# ------------------------------------------------------------

plot_df <- source_data %>%
  select(
    pair_id,
    region,
    bulk_V,
    log10_bulk_V,
    Bray_similarity,
    Jaccard_similarity
  ) %>%
  pivot_longer(
    cols = c(
      Bray_similarity,
      Jaccard_similarity
    ),
    names_to = "Metric",
    values_to = "Similarity"
  ) %>%
  filter(
    !is.na(bulk_V),
    bulk_V > 0,
    !is.na(Similarity)
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Bray_similarity" =
        "Bray-Curtis similarity",
      "Jaccard_similarity" =
        "Jaccard similarity"
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Bray-Curtis similarity",
        "Jaccard similarity"
      )
    )
  )


# ------------------------------------------------------------
# 14. Revised statistics labels
# ------------------------------------------------------------

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

format_r2 <- function(x) {
  if (is.na(x)) return("NA")
  sprintf("%.3f", x)
}

lab_bray <- sprintf(
  "Bray-Curtis: adj. R² = %s, P = %s",
  format_r2(res_bray$adjusted_R2),
  format_p(res_bray$permutation_P)
)

lab_jac <- sprintf(
  "Jaccard: adj. R² = %s, P = %s",
  format_r2(res_jac$adjusted_R2),
  format_p(res_jac$permutation_P)
)


# ------------------------------------------------------------
# 15. plot
# ------------------------------------------------------------

p_B <- ggplot(
  plot_df,
  aes(
    x = bulk_V,
    y = Similarity,
    color = Metric
  )
) +

  # points
  geom_point(
    size = 2.6,
    alpha = 0.74
  ) +

  # confidence ribbons
  geom_smooth(
    aes(fill = Metric),
    method = "lm",
    formula = y ~ log10(x),
    se = TRUE,
    linewidth = 0,
    alpha = 0.22,
    color = NA
  ) +

  # fitted lines
  geom_smooth(
    method = "lm",
    formula = y ~ log10(x),
    se = FALSE,
    linewidth = 1.15
  ) +

  scale_x_log10() +

  scale_color_manual(
    values = c(
      "Bray-Curtis similarity" = COL_BRAY_LINE,
      "Jaccard similarity" = COL_JAC_LINE
    )
  ) +

  scale_fill_manual(
    values = c(
      "Bray-Curtis similarity" = COL_BRAY_FILL,
      "Jaccard similarity" = COL_JAC_FILL
    )
  ) +

  labs(
    x = expression(
      "Bulk-soil V (mg kg"^{-1}*")"
    ),
    y = "Bulk-rhizosphere bacterial similarity",
    color = NULL,
    fill = NULL
  ) +

  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = lab_bray,
    hjust = -0.03,
    vjust = 1.35,
    size = 4.2,
    family = "Times New Roman",
    color = COL_BRAY_LINE
  ) +

  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = lab_jac,
    hjust = -0.03,
    vjust = 2.75,
    size = 4.2,
    family = "Times New Roman",
    color = COL_JAC_LINE
  ) +

  coord_cartesian(
    clip = "off"
  ) +

  theme_bw(
    base_family = "Times New Roman",
    base_size = 12.5
  ) +

  theme(
    panel.grid.major = element_line(
      color = "#D9D9D9",
      linewidth = 0.42
    ),
    panel.grid.minor = element_line(
      color = "#ECECEC",
      linewidth = 0.28
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.50
    ),
    axis.title = element_text(
      size = 13.5,
      color = "black"
    ),
    axis.text = element_text(
      size = 11.5,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.50,
      color = "black"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(
      size = 10.5
    ),
    legend.key.width = grid::unit(
      0.80,
      "cm"
    ),
    legend.margin = margin(
      t = -1,
      r = 0,
      b = 0,
      l = 0
    ),
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 6,
      l = 6
    )
  )


# ------------------------------------------------------------
# 16. Save plots
# ------------------------------------------------------------

ggsave(
  filename = OUT_PDF,
  plot = p_B,
  width = 5.0,
  height = 4.1,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = OUT_SVG,
  plot = p_B,
  width = 5.0,
  height = 4.1,
  units = "in",
  device = svglite::svglite
)

ggsave(
  filename = OUT_PNG,
  plot = p_B,
  width = 5.0,
  height = 4.1,
  units = "in",
  dpi = 600,
  bg = "white"
)

print(p_B)

cat("\n========================================\n")
cat("Fig. 4B finished.\n")
cat("All outputs saved to:\n", OUT_DIR, "\n")
cat("========================================\n")
