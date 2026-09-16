
# ============================================================
# FEAST pathway analysis across matched soil–plant systems
#
# Main paths:
#   1) Bulk soil -> Rhizosphere
#   2) Rhizosphere -> Endophyte
#   3) Endophyte -> Plant surface
#   4) Rhizosphere -> Plant surface
#   5) Bulk soil -> Plant surface
#
# ============================================================

script_dir <- function() {
  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      p <- frames[[i]]$ofile
      if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p)))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}

BASE_DIR <- script_dir()
LOCAL_LIBRARY <- file.path(dirname(BASE_DIR), "R_library")
if (dir.exists(LOCAL_LIBRARY)) .libPaths(c(LOCAL_LIBRARY, .libPaths()))
if (!requireNamespace("FEAST", quietly = TRUE)) {
  stop("FEAST is not available. Expected local library: ", LOCAL_LIBRARY, call. = FALSE)
}
setwd(BASE_DIR)

COUNT_FILE <- "FEAST_pathway_counts.tsv"
MAP_FILE   <- "FEAST_pathway_system_mapping.tsv"
OUTDIR     <- "FEAST_pathway_outputs"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

C <- read.table(
  COUNT_FILE, header = TRUE, row.names = 1, sep = "\t",
  check.names = FALSE, quote = "", comment.char = ""
)

map <- read.table(
  MAP_FILE, header = TRUE, sep = "\t",
  check.names = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE, na.strings = c("", "NA")
)

C <- as.matrix(C)
storage.mode(C) <- "numeric"

if (any(C < 0, na.rm = TRUE)) stop("Negative counts found.")
if (any(abs(C - round(C)) > 1e-8, na.rm = TRUE)) stop("Count matrix contains non-integer values.")
storage.mode(C) <- "integer"

# ---------- helpers ----------
nonempty <- function(x) {
  x <- x[!is.na(x)]
  x[nchar(x) > 0]
}

source_col_match <- function(raw, source_id) {
  # FEAST commonly returns "<SampleID>_<Env>"
  hits <- colnames(raw)[startsWith(colnames(raw), paste0(source_id, "_"))]
  if (length(hits) == 1) return(hits)
  hits2 <- colnames(raw)[grepl(source_id, colnames(raw), fixed = TRUE)]
  if (length(hits2) >= 1) return(hits2[1])
  return(NA_character_)
}

recover_sink_id <- function(row_name, sink_ids) {
  hit <- sink_ids[sapply(sink_ids, function(s) startsWith(row_name, paste0(s, "_")))]
  if (length(hit) == 1) return(hit)
  hit2 <- sink_ids[sapply(sink_ids, function(s) grepl(s, row_name, fixed = TRUE))]
  if (length(hit2) >= 1) return(hit2[1])
  return(row_name)
}

run_feast_model <- function(system_id, path_name, source_ids, sink_ids) {
  source_ids <- nonempty(source_ids)
  sink_ids   <- nonempty(sink_ids)

  if (length(source_ids) == 0 || length(sink_ids) == 0) {
    return(NULL)
  }

  missing_ids <- setdiff(c(source_ids, sink_ids), rownames(C))
  if (length(missing_ids) > 0) {
    stop(
      paste0(
        "Missing samples for System ", system_id, ", ", path_name, ": ",
        paste(missing_ids, collapse = ", ")
      )
    )
  }

  use_ids <- c(source_ids, sink_ids)
  Csub <- C[use_ids, , drop = FALSE]

  # Remove ASVs absent from this run only
  Csub <- Csub[, colSums(Csub) > 0, drop = FALSE]

  meta <- data.frame(
    Env = use_ids,
    SourceSink = c(rep("Source", length(source_ids)), rep("Sink", length(sink_ids))),
    id = c(rep(NA, length(source_ids)), seq_along(sink_ids)),
    row.names = use_ids,
    stringsAsFactors = FALSE
  )

  run_dir <- file.path(
    OUTDIR,
    paste0("System", system_id, "_", gsub("[^A-Za-z0-9]+", "_", path_name))
  )
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(run_dir)

  prefix <- "FEAST"

  FEAST::FEAST(
    C = Csub,
    metadata = meta,
    EM_iterations = 1000,
    COVERAGE = NULL,
    different_sources_flag = 0,
    dir_path = ".",
    outfile = prefix
  )

  raw_file <- paste0(prefix, "_source_contributions_matrix.txt")
  if (!file.exists(raw_file)) {
    stop(paste("Expected FEAST output not found:", file.path(run_dir, raw_file)))
  }

  raw <- read.table(
    raw_file, header = TRUE, row.names = 1, sep = "\t",
    check.names = FALSE, quote = "\"", comment.char = ""
  )

  # Match source columns
  source_cols <- setNames(
    vapply(source_ids, function(x) source_col_match(raw, x), character(1)),
    source_ids
  )

  if (any(is.na(source_cols))) {
    cat("\nAvailable FEAST columns:\n")
    print(colnames(raw))
    stop(
      paste(
        "Could not match source columns:",
        paste(names(source_cols)[is.na(source_cols)], collapse = ", ")
      )
    )
  }

  if (!("Unknown" %in% colnames(raw))) {
    stop("FEAST output has no 'Unknown' column.")
  }

  sink_recovered <- vapply(
    rownames(raw),
    recover_sink_id,
    character(1),
    sink_ids = sink_ids
  )

  # Total contribution from the source class specified for this path.
  # For Endophyte -> Plant surface, this sums ER/ES/EL source fractions.
  per_sink <- data.frame(
    System = system_id,
    Path = path_name,
    Sink = sink_recovered,
    SourceContribution = rowSums(raw[, source_cols, drop = FALSE], na.rm = TRUE),
    Unknown = as.numeric(raw[, "Unknown"]),
    stringsAsFactors = FALSE
  )

  per_sink$KnownOther <- pmax(
    0,
    1 - per_sink$SourceContribution - per_sink$Unknown
  )

  per_sink$SumCheck <- per_sink$SourceContribution + per_sink$Unknown + per_sink$KnownOther

  return(per_sink)
}

# ---------- run the five path models ----------
all_sink_results <- list()
k <- 1

for (i in seq_len(nrow(map))) {
  r <- map[i, ]
  sys <- r$System

  endo <- nonempty(c(r$Endophyte_root, r$Endophyte_stem, r$Endophyte_leaf))
  surf <- nonempty(c(r$Surface_root, r$Surface_stem, r$Surface_leaf))

  # 1) Bulk -> Rhizosphere
  z <- run_feast_model(
    sys, "Bulk_to_Rhizosphere",
    source_ids = r$Bulk,
    sink_ids = r$Rhizosphere
  )
  if (!is.null(z)) { all_sink_results[[k]] <- z; k <- k + 1 }

  # 2) Rhizosphere -> Endophyte
  z <- run_feast_model(
    sys, "Rhizosphere_to_Endophyte",
    source_ids = r$Rhizosphere,
    sink_ids = endo
  )
  if (!is.null(z)) { all_sink_results[[k]] <- z; k <- k + 1 }

  # 3) Endophyte -> Plant surface
  z <- run_feast_model(
    sys, "Endophyte_to_Plant_surface",
    source_ids = endo,
    sink_ids = surf
  )
  if (!is.null(z)) { all_sink_results[[k]] <- z; k <- k + 1 }

  # 4) Rhizosphere -> Plant surface
  z <- run_feast_model(
    sys, "Rhizosphere_to_Plant_surface",
    source_ids = r$Rhizosphere,
    sink_ids = surf
  )
  if (!is.null(z)) { all_sink_results[[k]] <- z; k <- k + 1 }

  # 5) Bulk -> Plant surface
  z <- run_feast_model(
    sys, "Bulk_to_Plant_surface",
    source_ids = r$Bulk,
    sink_ids = surf
  )
  if (!is.null(z)) { all_sink_results[[k]] <- z; k <- k + 1 }
}

per_sink <- do.call(rbind, all_sink_results)

# ---------- system-level means ----------
# Equal weighting of sinks within system for a given path.
per_system <- aggregate(
  cbind(SourceContribution, Unknown, KnownOther) ~ System + Path,
  data = per_sink,
  FUN = mean
)

# Add number of sinks used per system/path
n_sinks <- aggregate(
  Sink ~ System + Path,
  data = per_sink,
  FUN = length
)
colnames(n_sinks)[colnames(n_sinks) == "Sink"] <- "n_sinks"
per_system <- merge(per_system, n_sinks, by = c("System","Path"), all.x = TRUE)

# ---------- across-system summary ----------
path_levels <- c(
  "Bulk_to_Rhizosphere",
  "Rhizosphere_to_Endophyte",
  "Endophyte_to_Plant_surface",
  "Rhizosphere_to_Plant_surface",
  "Bulk_to_Plant_surface"
)

summary_list <- lapply(path_levels, function(pp) {
  x <- per_system$SourceContribution[per_system$Path == pp]
  x <- x[is.finite(x)]

  n <- length(x)
  mn <- if (n > 0) mean(x) else NA_real_
  sdv <- if (n > 1) sd(x) else NA_real_
  se <- if (n > 1) sdv / sqrt(n) else NA_real_
  med <- if (n > 0) median(x) else NA_real_
  ci_lo <- if (n > 1) mn - qt(0.975, df = n - 1) * se else NA_real_
  ci_hi <- if (n > 1) mn + qt(0.975, df = n - 1) * se else NA_real_

  data.frame(
    Path = pp,
    n_systems = n,
    MeanContribution = mn,
    SD = sdv,
    SE = se,
    Median = med,
    CI95_low = ci_lo,
    CI95_high = ci_hi,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_list)

# Friendly path labels for plotting/diagram
label_map <- c(
  Bulk_to_Rhizosphere = "Bulk soil → Rhizosphere",
  Rhizosphere_to_Endophyte = "Rhizosphere → Endophyte",
  Endophyte_to_Plant_surface = "Endophyte → Plant surface",
  Rhizosphere_to_Plant_surface = "Rhizosphere → Plant surface",
  Bulk_to_Plant_surface = "Bulk soil → Plant surface"
)

summary_df$PathLabel <- unname(label_map[summary_df$Path])
summary_df$MeanPercent <- 100 * summary_df$MeanContribution
summary_df$SEPercent <- 100 * summary_df$SE
summary_df$CI95_low_percent <- 100 * summary_df$CI95_low
summary_df$CI95_high_percent <- 100 * summary_df$CI95_high

# Reorder columns
summary_df <- summary_df[, c(
  "Path", "PathLabel", "n_systems",
  "MeanContribution", "MeanPercent",
  "SD", "SE", "SEPercent",
  "Median", "CI95_low", "CI95_high",
  "CI95_low_percent", "CI95_high_percent"
)]

# ---------- outputs ----------
write.table(
  per_sink,
  file.path(OUTDIR, "FEAST_pathway_per_sink.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  per_system,
  file.path(OUTDIR, "FEAST_pathway_per_system.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  summary_df,
  file.path(OUTDIR, "FEAST_pathway_mean_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.csv(
  summary_df,
  file.path(OUTDIR, "FEAST_pathway_mean_summary.csv"),
  row.names = FALSE
)

# A compact label file for directly annotating the pathway diagram
diagram_labels <- summary_df[, c(
  "PathLabel", "n_systems", "MeanPercent", "SEPercent",
  "CI95_low_percent", "CI95_high_percent"
)]

write.csv(
  diagram_labels,
  file.path(OUTDIR, "FEAST_pathway_diagram_labels.csv"),
  row.names = FALSE
)

cat("\n========================================\n")
cat("FEAST pathway analysis completed.\n")
cat("Main summary (equal system weight):\n")
print(summary_df[, c("PathLabel","n_systems","MeanPercent","SEPercent")], row.names = FALSE)
cat("\nKey file for the diagram:\n")
cat(file.path(OUTDIR, "FEAST_pathway_diagram_labels.csv"), "\n")
cat("========================================\n")
