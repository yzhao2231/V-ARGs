# ============================================================
# Fig. 4C — Bulk-soil V vs soil-derived microbial contribution
# ALL PLANT SINKS version
# ============================================================
#
# Required raw files in this analysis folder (or subfolders):
#   Soil geochemistry.xlsx
#   FEAST_v2_full_14source_contributions.csv
#
# All plant sinks are included:
#   SR, SS, SL, ER, ES, EL (where available)
#
# Soil-derived microbial contribution for each sink =
#   sum of all 14 identified soil-source FEAST weights
#   = 1 - Unknown
#
# Within each matched system, all available plant sinks are
# equally averaged to obtain one system-level contribution.
#
# Formal test:
#   exact Spearman rank correlation
#   + exact within-region permutation sensitivity analysis
#
# Visual line:
#   simple linear fit + 95% CI (visualization only)
#
# Outputs -> ./output/
# ============================================================

required_pkgs <- c("readxl","dplyr","ggplot2","svglite","rstudioapi")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)
]
if(length(missing_pkgs)>0){
  stop("Install missing packages: ",
       paste(missing_pkgs, collapse=", "))
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(svglite)
})

get_script_dir <- function(){
  if(requireNamespace("rstudioapi", quietly=TRUE) &&
     rstudioapi::isAvailable()){
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error=function(e) "")
    if(nzchar(p)){
      return(dirname(normalizePath(p, winslash="/", mustWork=TRUE)))
    }
  }
  normalizePath(getwd(), winslash="/", mustWork=TRUE)
}

BASE_DIR <- get_script_dir()
OUT_DIR <- file.path(BASE_DIR,"output")
dir.create(OUT_DIR,showWarnings=FALSE,recursive=TRUE)

geo_hits <- list.files(
  BASE_DIR,
  pattern="^Soil geochemistry\\.xlsx$",
  recursive=TRUE, full.names=TRUE, ignore.case=TRUE
)
feast_hits <- list.files(
  BASE_DIR,
  pattern="^FEAST_v2_full_14source_contributions\\.csv$",
  recursive=TRUE, full.names=TRUE, ignore.case=TRUE
)

if(length(geo_hits)==0) stop("Cannot find Soil geochemistry.xlsx")
if(length(feast_hits)==0) stop("Cannot find FEAST_v2_full_14source_contributions.csv")

GEOCHEM_FILE <- geo_hits[1]
FEAST_FILE <- feast_hits[1]

system_map <- data.frame(
  system=1:7,
  region=c("Chengde","Chengde","Huaihua","Huaihua","Huaihua","Panzhihua","Panzhihua"),
  crop=c("Corn","Sweet potato","Corn","Sweet potato","Peanut","Corn","Sweet potato"),
  bulk_sample=c("A2-11","A2-10","B2-8","B2-1","B2-11","C2-17","C2-13"),
  stringsAsFactors=FALSE
)

# -------------------------
# Read V
# -------------------------
geo <- as.data.frame(read_excel(GEOCHEM_FILE, sheet=1, .name_repair="minimal"))
names(geo) <- trimws(names(geo))

if(!all(c("ID","V") %in% names(geo))){
  stop("Geochemistry workbook must contain ID and V columns.")
}

geo_v <- geo %>%
  transmute(
    sample_id=trimws(as.character(ID)),
    bulk_soil_V_mg_kg=suppressWarnings(as.numeric(V))
  ) %>%
  distinct(sample_id,.keep_all=TRUE)

# -------------------------
# Read microbial FEAST
# -------------------------
feast <- read.csv(
  FEAST_FILE,
  check.names=FALSE,
  stringsAsFactors=FALSE
)
names(feast) <- trimws(names(feast))

if(!all(c("Sink","Unknown") %in% names(feast))){
  stop("FEAST file must contain Sink and Unknown columns.")
}

source_cols <- setdiff(names(feast),c("Sink","Unknown"))
if(length(source_cols)!=14){
  stop("Expected 14 soil-source columns, detected ",length(source_cols))
}

for(nm in c(source_cols,"Unknown")){
  feast[[nm]] <- suppressWarnings(as.numeric(feast[[nm]]))
}

# standardize FEAST scale to fraction
row_totals <- rowSums(feast[,c(source_cols,"Unknown"),drop=FALSE],na.rm=TRUE)
if(median(row_totals,na.rm=TRUE)>1.5){
  feast[,c(source_cols,"Unknown")] <-
    feast[,c(source_cols,"Unknown")] / 100
}

# parse system from SR-1, SS-1, ...
feast <- feast %>%
  mutate(
    system=as.integer(sub(".*-","",Sink)),
    sink_type=sub("-.*","",Sink),
    soil_derived_fraction=
      rowSums(across(all_of(source_cols)),na.rm=TRUE),
    soil_derived_percent=soil_derived_fraction*100,
    unknown_percent=Unknown*100
  )

# -------------------------
# ALL available plant sinks
# -------------------------
all_plant_sink <- feast %>%
  filter(sink_type %in% c("SR","SS","SL","ER","ES","EL"))

system_contribution <- all_plant_sink %>%
  group_by(system) %>%
  summarise(
    n_plant_sinks=n(),
    plant_sinks=paste(Sink,collapse="; "),
    soil_derived_microbe_contribution_percent=
      mean(soil_derived_percent,na.rm=TRUE),
    mean_unknown_percent=
      mean(unknown_percent,na.rm=TRUE),
    .groups="drop"
  )

system_data <- system_map %>%
  left_join(system_contribution,by="system") %>%
  left_join(
    geo_v %>% rename(bulk_sample=sample_id),
    by="bulk_sample"
  )

write.csv(
  all_plant_sink,
  file.path(OUT_DIR,"Fig4C_all_plant_sinks_sink_level.csv"),
  row.names=FALSE,fileEncoding="UTF-8"
)

write.csv(
  system_data,
  file.path(OUT_DIR,"Fig4C_all_plant_sinks_source_data.csv"),
  row.names=FALSE,fileEncoding="UTF-8"
)

cat("\nSystem-level data:\n")
print(system_data)

# -------------------------
# Exact permutation helpers
# -------------------------
all_permutations <- function(x){
  if(length(x)==1) return(matrix(x,nrow=1))
  do.call(rbind,lapply(seq_along(x),function(i){
    cbind(x[i],all_permutations(x[-i]))
  }))
}

rho_fun <- function(x,y){
  suppressWarnings(cor(x,y,method="spearman",use="complete.obs"))
}

x <- system_data$bulk_soil_V_mg_kg
y <- system_data$soil_derived_microbe_contribution_percent

rho_obs <- rho_fun(x,y)

# Global exact 7! permutations
perm_global <- all_permutations(x)
rho_global <- apply(
  perm_global,1,function(xx) rho_fun(xx,y)
)
p_global <- mean(
  abs(rho_global) >= abs(rho_obs)-1e-12
)

# Within-region exact permutations: 2! x 3! x 2! = 24
region_idx <- split(seq_len(nrow(system_data)),system_data$region)
region_perm <- lapply(region_idx,function(idx) all_permutations(x[idx]))

perm_grid <- expand.grid(
  lapply(region_perm,function(m) seq_len(nrow(m)))
)

rho_region <- numeric(nrow(perm_grid))

for(i in seq_len(nrow(perm_grid))){
  xp <- x
  for(g in seq_along(region_idx)){
    idx <- region_idx[[g]]
    xp[idx] <- region_perm[[g]][perm_grid[i,g],,drop=TRUE]
  }
  rho_region[i] <- rho_fun(xp,y)
}

p_region <- mean(
  abs(rho_region) >= abs(rho_obs)-1e-12
)

# Visual LM only
visual_lm <- lm(y~x)
visual_sm <- summary(visual_lm)

stats_out <- data.frame(
  N_systems=length(x),
  Spearman_rho=rho_obs,
  Global_exact_P=p_global,
  Within_region_exact_P=p_region,
  Global_permutations=nrow(perm_global),
  Within_region_permutations=nrow(perm_grid),
  Visual_LM_R2=visual_sm$r.squared,
  Visual_LM_P=coef(visual_sm)[2,"Pr(>|t|)"]
)

write.csv(
  stats_out,
  file.path(OUT_DIR,"Fig4C_all_plant_sinks_statistics.csv"),
  row.names=FALSE
)

cat("\nStatistics:\n")
print(stats_out)

# -------------------------
# Plot
# -------------------------
COL_POINT <- "#B67F52"
COL_LINE  <- "#9A6844"
COL_CI    <- "#DBB999"

format_p <- function(p){
  if(p<0.001) "< 0.001" else sprintf("%.4f",p)
}

stat_label <- sprintf(
  "\u03c1 = %.3f, exact P = %s",
  rho_obs,
  format_p(p_global)
)

p_C <- ggplot(
  system_data,
  aes(
    x=bulk_soil_V_mg_kg,
    y=soil_derived_microbe_contribution_percent
  )
) +
  geom_smooth(
    method="lm",
    formula=y~x,
    se=TRUE,
    color=COL_LINE,
    fill=COL_CI,
    linewidth=1.15,
    alpha=0.32
  ) +
  geom_point(
    shape=21,
    size=4.0,
    stroke=0.8,
    color=COL_LINE,
    fill=COL_POINT
  ) +
  labs(
    x=expression("Bulk-soil V (mg kg"^{-1}*")"),
    y="Soil-derived microbial contribution (%)"
  ) +
  annotate(
    "text",
    x=-Inf,y=Inf,
    label=stat_label,
    hjust=-0.05,vjust=1.30,
    family="Arial",
    size=4.8,
    color="black"
  ) +
  coord_cartesian(clip="off") +
  theme_bw(base_family="Times New Roman",base_size=13) +
  theme(
    panel.grid.major=element_line(color="#D9D9D9",linewidth=0.42),
    panel.grid.minor=element_line(color="#ECECEC",linewidth=0.28),
    panel.border=element_rect(color="black",fill=NA,linewidth=0.80),
    axis.title=element_text(size=14,color="black"),
    axis.text=element_text(size=12,color="black"),
    axis.ticks=element_line(color="black",linewidth=0.50),
    plot.margin=margin(8,8,6,6)
  )

print(p_C)

ggsave(
  file.path(OUT_DIR,"Fig4C_all_plant_sinks.pdf"),
  p_C,width=4.6,height=4.0,units="in",
  device=cairo_pdf
)

ggsave(
  file.path(OUT_DIR,"Fig4C_all_plant_sinks.svg"),
  p_C,width=4.6,height=4.0,units="in",
  device=svglite::svglite
)

ggsave(
  file.path(OUT_DIR,"Fig4C_all_plant_sinks.png"),
  p_C,width=4.6,height=4.0,units="in",
  dpi=600,bg="white"
)

cat("\nDone. Outputs saved to:\n",OUT_DIR,"\n")
