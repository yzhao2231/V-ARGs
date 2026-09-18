# ============================================================
# Fig. 4D — Microbial connectivity vs soil–plant resistome connectivity
# Formal model
# ============================================================
#
# Total ARG matched connectivity index (MCI):
#   MCI = S(P, matched-soil centroid) - S(P, unmatched-soil centroid)
#   where S = 1 - Bray-Curtis dissimilarity.
#
# Formal model:
#   MCI_z ~ microbial_contribution_z + System + Habitat

# Significance:
#   Freedman-Lane residual permutation; residuals permuted within System.
#   B = 99,999 permutations.
#
# Plot:
#   Raw sink-level points + simple lm fit and 95% CI for visualization.
#   Formal inference is from the adjusted permutation model above.
#
# ============================================================

required_pkgs <- c("ggplot2", "svglite")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

# ------------------------------------------------------------
# 1. Locate analysis folder
# ------------------------------------------------------------
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    if (nzchar(p)) {
      return(dirname(normalizePath(p, winslash = "/", mustWork = TRUE)))
    }
  }

  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    p <- sub("^--file=", "", hit[1])
    return(dirname(normalizePath(p, winslash = "/", mustWork = TRUE)))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

BASE_DIR <- get_script_dir()
OUT_DIR <- file.path(BASE_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

find_one <- function(filename) {
  hits <- list.files(
    BASE_DIR,
    pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) == 0) stop("Cannot find: ", filename)
  hits[1]
}

ARG_FILE   <- find_one("01 ARG_MGE_relative_abundance.csv")
FEAST_FILE <- find_one("FEAST_v2_full_14source_contributions.csv")

# ------------------------------------------------------------
# 2. Matched soil–plant systems
# ------------------------------------------------------------
system_map <- data.frame(
  System = 1:7,
  Region = c(
    "Chengde", "Chengde",
    "Huaihua", "Huaihua", "Huaihua",
    "Panzhihua", "Panzhihua"
  ),
  Rhizosphere_sample = c(
    "A1-5", "A1-12", "B1-13", "B1-9",
    "B1-7", "C1-11", "C1-16"
  ),
  Bulk_soil_sample = c(
    "A2-11", "A2-10", "B2-8", "B2-1",
    "B2-11", "C2-17", "C2-13"
  ),
  stringsAsFactors = FALSE
)

all_soils <- c(
  system_map$Rhizosphere_sample,
  system_map$Bulk_soil_sample
)

# ------------------------------------------------------------
# 3. Read ARG relative-abundance matrix
# ------------------------------------------------------------
arg <- read.csv(
  ARG_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(c("Gene_type", "Assay") %in% names(arg))) {
  stop("ARG file must contain Gene_type and Assay columns.")
}

arg <- arg[trimws(as.character(arg$Gene_type)) == "ARG", , drop = FALSE]

if (nrow(arg) == 0) stop("No ARG rows detected.")

get_profile <- function(sample_id) {
  if (!sample_id %in% names(arg)) {
    stop("Sample not found in ARG matrix: ", sample_id)
  }
  x <- suppressWarnings(as.numeric(arg[[sample_id]]))
  x[is.na(x)] <- 0
  x
}

get_centroid <- function(sample_ids) {
  miss <- setdiff(sample_ids, names(arg))
  if (length(miss) > 0) {
    stop("Samples missing from ARG matrix: ", paste(miss, collapse = ", "))
  }

  x <- sapply(sample_ids, get_profile)

  if (is.null(dim(x))) {
    return(as.numeric(x))
  }

  rowMeans(x, na.rm = TRUE)
}

bray_similarity <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)

  denominator <- sum(a + b, na.rm = TRUE)

  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }

  1 - sum(abs(a - b), na.rm = TRUE) / denominator
}

# ------------------------------------------------------------
# 4. Read 14-source microbial FEAST output
# ------------------------------------------------------------
feast <- read.csv(
  FEAST_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(c("Sink", "Unknown") %in% names(feast))) {
  stop("FEAST file must contain Sink and Unknown columns.")
}

source_cols <- setdiff(names(feast), c("Sink", "Unknown"))

if (length(source_cols) != 14) {
  stop("Expected 14 identified soil-source columns; detected ",
       length(source_cols), ".")
}

for (nm in c(source_cols, "Unknown")) {
  feast[[nm]] <- suppressWarnings(as.numeric(feast[[nm]]))
}

row_totals <- rowSums(
  feast[, c(source_cols, "Unknown"), drop = FALSE],
  na.rm = TRUE
)

if (median(row_totals, na.rm = TRUE) > 1.5) {
  feast[, c(source_cols, "Unknown")] <-
    feast[, c(source_cols, "Unknown")] / 100
}

# ------------------------------------------------------------
# 5. Formal model
#    Tissue is intentionally NOT included.
# ------------------------------------------------------------
z <- function(x) as.numeric(scale(x))

d <- source_data
d$MCI_z <- z(d$Total_ARG_MCI)
d$Microbial_contribution_z <-
  z(d$Soil_derived_microbial_contribution_percent)

d$System <- factor(d$System)
d$Habitat <- factor(d$Habitat, levels = c("Surface", "Endophyte"))

full_formula <-
  MCI_z ~ Microbial_contribution_z + System + Habitat

reduced_formula <-
  MCI_z ~ System + Habitat

full_model <- lm(full_formula, data = d)
reduced_model <- lm(reduced_formula, data = d)

beta_obs <- unname(
  coef(full_model)["Microbial_contribution_z"]
)

# ------------------------------------------------------------
# 6. Freedman-Lane residual permutation
#    Residuals permuted only within matched System.
# ------------------------------------------------------------
B <- 99999
set.seed(20260918)

X_full <- model.matrix(full_model)
X_red  <- model.matrix(reduced_model)

y <- d$MCI_z

coef_map <- solve(crossprod(X_full), t(X_full))

beta_row <- which(
  colnames(X_full) == "Microbial_contribution_z"
)

fitted_red <- as.numeric(
  X_red %*% coef(reduced_model)
)
resid_red <- resid(reduced_model)

system_index <- split(seq_len(nrow(d)), d$System)

beta_perm <- numeric(B)

for (b in seq_len(B)) {

  perm_resid <- resid_red

  for (idx in system_index) {
    perm_resid[idx] <- sample(resid_red[idx], replace = FALSE)
  }

  y_perm <- fitted_red + perm_resid

  beta_perm[b] <- sum(
    coef_map[beta_row, ] * y_perm
  )
}

perm_p <- (
  sum(abs(beta_perm) >= abs(beta_obs) - 1e-12) + 1
) / (B + 1)

# ------------------------------------------------------------
# 7. Supporting descriptive statistics
# ------------------------------------------------------------
spearman_test <- suppressWarnings(
  cor.test(
    source_data$Soil_derived_microbial_contribution_percent,
    source_data$Total_ARG_MCI,
    method = "spearman",
    exact = FALSE
  )
)

x_resid <- resid(
  lm(
    Microbial_contribution_z ~ System + Habitat,
    data = d
  )
)

y_resid <- resid(
  lm(
    MCI_z ~ System + Habitat,
    data = d
  )
)

partial_r <- cor(x_resid, y_resid, method = "pearson")

stats_out <- data.frame(
  N_plant_sinks = nrow(d),
  Standardized_beta = beta_obs,
  Freedman_Lane_permutations = B,
  Within_system_permutation_P = perm_p,
  Partial_Pearson_r = partial_r,
  Raw_Spearman_rho = unname(spearman_test$estimate),
  Raw_Spearman_asymptotic_P = spearman_test$p.value,
  Adjusted_R2_full_model = summary(full_model)$adj.r.squared,
  stringsAsFactors = FALSE
)

write.csv(
  stats_out,
  file.path(OUT_DIR, "Fig4D_statistics_system_habitat.csv"),
  row.names = FALSE
)

cat("\nFig. 4D formal statistics:\n")
print(stats_out)

# ------------------------------------------------------------
# 9. Plot
# ------------------------------------------------------------
COL_POINT <- "#B27697"
COL_LINE  <- "#8F5E7B"
COL_CI    <- "#D9B6C9"

stat_label <- sprintf(
  "\u03b2 = %.2f, permutation P %s",
  beta_obs,
  ifelse(
    perm_p < 0.001,
    "< 0.001",
    paste0("= ", sprintf("%.3f", perm_p))
  )
)

p_D <- ggplot(
  source_data,
  aes(
    x = Soil_derived_microbial_contribution_percent,
    y = Total_ARG_MCI
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = COL_LINE,
    fill = COL_CI,
    linewidth = 1.15,
    alpha = 0.30
  ) +
  geom_point(
    shape = 21,
    size = 3.5,
    stroke = 0.75,
    color = COL_LINE,
    fill = COL_POINT,
    alpha = 0.88
  ) +
  labs(
    x = "Soil-derived microbial contribution (%)",
    y = "Total ARG matched connectivity index (MCI)"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = stat_label,
    hjust = -0.05,
    vjust = 1.30,
    family = "Arial",
    size = 4.8,
    color = "black"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(
    base_family = "Arial",
    base_size = 13
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
      linewidth = 0.80
    ),
    axis.title = element_text(
      size = 14,
      color = "black"
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.50
    ),
    plot.margin = margin(8, 8, 6, 6)
  )

print(p_D)

ggsave(
  file.path(OUT_DIR, "Fig4D_plot_system_habitat.pdf"),
  p_D,
  width = 4.6,
  height = 4.0,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(OUT_DIR, "Fig4D_plot_system_habitat.svg"),
  p_D,
  width = 4.6,
  height = 4.0,
  units = "in",
  device = svglite::svglite
)

ggsave(
  file.path(OUT_DIR, "Fig4D_plot_system_habitat.png"),
  p_D,
  width = 4.6,
  height = 4.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

cat("\nDone. Outputs saved to:\n", OUT_DIR, "\n")
