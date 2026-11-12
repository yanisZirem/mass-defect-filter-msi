#!/usr/bin/env Rscript
# =============================================================================
# Mass defect filter (MDF) for MALDI-MSI ion lists
# Separates lipid-like ions from non-lipid ions (matrix, salts, peptides, ...)
# based on the mass defect of theoretical lipid ions.
#
# Usage:
#   Rscript mass_defect_filter.R ion_list.xlsx [data.imzML]
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(writexl)
})

# ----------------------------- CONFIG ----------------------------------------
args <- commandArgs(trailingOnly = TRUE)

xlsx_file <- if (length(args) >= 1) args[1] else "ion_list.xlsx"
imzml_file <- if (length(args) >= 2) args[2] else NA_character_  # optional
out_dir   <- "mass_defect_output"

col_mz         <- 1   # column holding the m/z values
col_annotation <- 7   # "Standardized annotation"
col_msms       <- 8   # "MS/MS confirmed (this study)"
first_row      <- 3   # first data row (row 1 = title, row 2 = headers)

md_tol   <- 0.010   # filter tolerance (Da) around the lipid envelope
ppm_tol  <- 10      # peak-matching tolerance (ppm) in the imzML
snr_min  <- 3       # S/N threshold for peak detection (imzML part)
mz_range <- c(400, 1000)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- MASSES ----------------------------------------
mono <- c(C = 12, H = 1.00782503207, N = 14.0030740048,
          O = 15.99491461956, P = 30.97376163, S = 31.97207100)
M_CH2  <- 12 + 2 * mono[["H"]]            # 14.01565
M_H2   <- 2 * mono[["H"]]                 # 2.01565
M_H2O  <- 2 * mono[["H"]] + mono[["O"]]
M_O    <- mono[["O"]]
e_mass <- 0.00054858                      # electron mass
M_Hp   <- mono[["H"]] - e_mass            # proton
M_Nap  <- 22.98976928 - e_mass
M_Kp   <- 38.96370668 - e_mass
M_NH4p <- mono[["N"]] + 4 * mono[["H"]] - e_mass

parse_formula <- function(s) {
  tok <- regmatches(s, gregexpr("[A-Z][a-z]?[0-9]*", s))[[1]]
  el  <- sub("[0-9]*$", "", tok)
  n   <- suppressWarnings(as.numeric(sub("^[A-Z][a-z]?", "", tok)))
  n[is.na(n)] <- 1
  sum(mono[el] * n)
}

# ---------------------- THEORETICAL LIPIDS -----------------------------------
# For each class: formula of a reference species (Cref, DBref);
# M(C, DB, nO) = Mref + (C-Cref)*CH2 - (DB-DBref)*H2 + nO*O
# (C, DB = total chain carbons and double bonds, LIPID MAPS convention)
classes <- data.frame(
  class = c("FA","LPC","LPE","LPS","LPI","LPA","LPG",
            "PC","PE","PS","PI","PA","PG","PC-O","PE-O",
            "SM","Cer","HexCer","DG","TG","CE","ST","ST-sulfate"),
  formula = c("C16H32O2","C24H50NO7P","C21H44NO7P","C22H44NO9P","C25H49O12P","C19H39O7P","C22H45O9P",
              "C40H80NO8P","C37H74NO8P","C38H74NO10P","C41H79O13P","C35H69O8P","C38H75O10P","C40H82NO7P","C37H76NO7P",
              "C39H79N2O6P","C34H67NO3","C40H77NO8","C35H68O5","C51H98O6","C45H78O2","C27H46O","C27H46O4S"),
  Cref  = c(16,16,16,16,16,16,16, 32,32,32,32,32,32,32,32, 34,34,34, 32,48,18, 27,27),
  DBref = c( 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1,  0, 0, 1,  0, 0),
  Cmin  = c(12,12,12,12,12,12,12, 26,26,26,26,26,26,26,26, 30,30,30, 28,42,12, 18,18),
  Cmax  = c(26,26,26,26,26,26,26, 46,46,46,46,46,46,46,46, 44,44,44, 42,62,24, 30,30),
  DBmax = c( 6, 6, 6, 6, 6, 6, 6, 12,12,12,12,12,12,12,12,  4, 4, 4,  8,12, 6,  6, 6),
  nOmax = c( 2, 1, 1, 1, 1, 1, 1,  2, 2, 2, 2, 2, 2, 2, 2,  2, 2, 2,  2, 2, 1,  5, 3),
  stringsAsFactors = FALSE
)
classes$Mref <- vapply(classes$formula, parse_formula, numeric(1))

adducts <- c("[M+H]+"      = M_Hp,
             "[M+Na]+"     = M_Nap,
             "[M+K]+"      = M_Kp,
             "[M+NH4]+"    = M_NH4p,
             "[M+H-H2O]+"  = M_Hp - M_H2O,
             "[M+2Na-H]+"  = 2 * M_Nap - M_Hp)

build_theoretical <- function() {
  out <- list()
  for (i in seq_len(nrow(classes))) {
    k <- classes[i, ]
    g <- expand.grid(C = k$Cmin:k$Cmax, DB = 0:k$DBmax, nO = 0:k$nOmax)
    g <- g[g$DB <= g$C / 2, ]                       # chemical sanity check
    M <- k$Mref + (g$C - k$Cref) * M_CH2 - (g$DB - k$DBref) * M_H2 + g$nO * M_O
    for (a in names(adducts)) {
      out[[length(out) + 1]] <- data.frame(class = k$class, mz = M + adducts[[a]])
    }
  }
  th <- do.call(rbind, out)
  th[th$mz >= mz_range[1] - 5 & th$mz <= mz_range[2] + 5, ]
}
theo <- build_theoretical()
theo$N   <- floor(theo$mz)          # nominal mass (floor convention)
theo$md  <- theo$mz - theo$N        # mass defect (0 <= md < 1)

# Lipid envelope per nominal mass (min/max mass defect)
env <- aggregate(md ~ N, data = theo, FUN = function(x) c(lo = min(x), hi = max(x)))
env <- data.frame(N = env$N, lo = env$md[, "lo"], hi = env$md[, "hi"])
# Nominal masses without any theoretical lipid reuse the envelope of neighbouring masses
all_N <- seq(floor(mz_range[1]) - 2, ceiling(mz_range[2]) + 2)
env <- merge(data.frame(N = all_N), env, by = "N", all.x = TRUE)
for (j in which(is.na(env$lo))) {
  nb <- env[abs(env$N - env$N[j]) <= 3 & !is.na(env$lo), ]
  env$lo[j] <- min(nb$lo); env$hi[j] <- max(nb$hi)
}

# ----------------------------- CLASSIFICATION --------------------------------
mdf_classify <- function(mz, tol = md_tol) {
  N  <- floor(mz)
  md <- mz - N
  lo <- env$lo[match(N, env$N)] - tol
  hi <- env$hi[match(N, env$N)] + tol
  cat <- ifelse(md < lo, "Non-lipid (defect below lipid window)",
         ifelse(md > hi, "Non-lipid (defect above lipid window)", "Lipid-like"))
  data.frame(mz = mz, nominal_mass = N, mass_defect = md,
             rel_mass_defect_ppm = md / mz * 1e6,
             window_lo = lo, window_hi = hi, mdf_class = cat,
             stringsAsFactors = FALSE)
}

# Kendrick (CH2 base): lipids of the same class share a constant KMD
kendrick <- function(mz) {
  km <- mz * 14 / M_CH2
  data.frame(kendrick_mass = km, KMD = round(km) - km)
}

# Clusters de matrice DHB (C7H6O4) : [nDHB - m H2O + adduit]+  (indicative)
dhb <- parse_formula("C7H6O4")
dhb_clusters <- function(tol_da = 0.005) {
  g <- expand.grid(n = 1:6, m = 0:6)
  g <- g[g$m <= g$n, ]
  out <- list()
  for (a in c("[M+H]+", "[M+Na]+", "[M+K]+")) {
    out[[a]] <- data.frame(label = sprintf("[%dDHB-%dH2O%s", g$n, g$m, sub("\\[M", "", a)),
                           mz = g$n * dhb - g$m * M_H2O + adducts[[a]])
  }
  do.call(rbind, out)
}
dhb_tab <- dhb_clusters()
match_dhb <- function(mz, tol_da = 0.005) {
  vapply(mz, function(x) {
    i <- which(abs(dhb_tab$mz - x) <= tol_da)
    if (length(i)) paste(dhb_tab$label[i], collapse = " / ") else NA_character_
  }, character(1))
}

# ----------------------------- LECTURE EXCEL ---------------------------------
raw <- read_excel(xlsx_file, col_names = FALSE, skip = first_row - 1)
d <- data.frame(mz = suppressWarnings(as.numeric(raw[[col_mz]])),
                annotation = as.character(raw[[col_annotation]]),
                msms_confirmed = as.character(raw[[col_msms]]),
                stringsAsFactors = FALSE)
d <- d[!is.na(d$mz), ]          # drops the "Total" row
d$annotated <- !is.na(d$annotation)

res <- cbind(d, mdf_classify(d$mz)[, -1], kendrick(d$mz))
res$possible_DHB_cluster <- match_dhb(d$mz)
res$lipid_like <- res$mdf_class == "Lipid-like"

# Consistency between annotations and filter
res$check <- ifelse(res$annotated & !res$lipid_like, "Annotated but fails MDF (to re-check)",
             ifelse(!res$annotated & res$lipid_like, "Unannotated, lipid-like (candidate)",
             ifelse(res$annotated, "Annotated, passes MDF", "Unannotated, non-lipid")))

# ----------------------------- SUMMARY ---------------------------------------
summary_tab <- as.data.frame(table(res$mdf_class), stringsAsFactors = FALSE)
names(summary_tab) <- c("Category", "n_ions")
summary_tab$percent <- round(100 * summary_tab$n_ions / sum(summary_tab$n_ions), 1)
check_tab <- as.data.frame(table(res$check), stringsAsFactors = FALSE)
names(check_tab) <- c("Annotation vs MDF", "n_ions")

cat("\n==== Mass defect filter (tolerance +/-", md_tol, "Da) ====\n")
print(summary_tab, row.names = FALSE)
cat("\n")
print(check_tab, row.names = FALSE)
cat("\nAnnotated ions failing the filter:\n")
print(res[res$check == "Annotated but fails MDF (to re-check)",
          c("mz", "annotation", "mass_defect", "window_lo", "window_hi")], row.names = FALSE)

# ----------------------------- FIGURES ---------------------------------------
env_plot <- env[env$N >= mz_range[1] - 1 & env$N <= mz_range[2], ]
env_plot <- data.frame(mz = c(env_plot$N, env_plot$N + 1),
                       lo = rep(env_plot$lo - md_tol, 2),
                       hi = rep(env_plot$hi + md_tol, 2))
env_plot <- env_plot[order(env_plot$mz), ]

cols <- c("Lipid-like" = "#1b7837",
          "Non-lipid (defect below lipid window)" = "#c51b7d",
          "Non-lipid (defect above lipid window)" = "#e08214")

p1 <- ggplot() +
  geom_ribbon(data = env_plot, aes(mz, ymin = lo, ymax = hi),
              fill = "#1b7837", alpha = 0.12) +
  geom_point(data = res, aes(mz, mass_defect, colour = mdf_class,
                             shape = annotated), size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                     labels = c("unannotated", "annotated"), name = NULL) +
  labs(x = expression(italic(m/z)), y = "Mass defect (Da)",
       title = "Mass defect filter - DHB, positive mode",
       subtitle = "Green band: theoretical lipid envelope (FA, GPL, SM, Cer, DG, TG, CE, sterols; H/Na/K/NH4 adducts)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "mass_defect_plot.png"), p1, width = 8.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "mass_defect_plot.pdf"), p1, width = 8.5, height = 5.5)

p2 <- ggplot(res, aes(kendrick_mass, KMD, colour = mdf_class, shape = annotated)) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                     labels = c("unannotated", "annotated"), name = NULL) +
  labs(x = "Kendrick mass (CH2)", y = "Kendrick mass defect",
       title = "Kendrick mass defect plot") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "kendrick_plot.png"), p2, width = 8.5, height = 5.5, dpi = 300)

# ----------------------------- EXPORT ----------------------------------------
write_xlsx(list(MDF_all_ions = res,
                Non_lipid_ions = res[!res$lipid_like, ],
                Summary = summary_tab,
                Annotation_vs_MDF = check_tab,
                Lipid_envelope = env),
           file.path(out_dir, "mass_defect_results.xlsx"))

# ============================ imzML PART (optional) ==========================
# Applies the filter (1) to the listed ions and (2) to ALL peaks detected in the
# mean spectrum, and reports the lipid-like vs non-lipid intensity share.
if (!is.na(imzml_file) && file.exists(imzml_file)) {
  suppressPackageStartupMessages(library(Cardinal))
  cat("\nReading", imzml_file, "...\n")
  msi <- if (exists("readMSIData", asNamespace("Cardinal"))) readMSIData(imzml_file) else readImzML(imzml_file)

  mzs <- as.numeric(Cardinal::mz(msi))
  # mean spectrum (Cardinal 3.x; matrix fallback otherwise)
  mean_int <- tryCatch({
    s <- summarizeFeatures(msi, "mean")
    as.numeric(if (is.data.frame(s) || inherits(s, "DataFrame")) s$mean else fData(s)$mean)
  }, error = function(e) rowMeans(as.matrix(Cardinal::spectra(msi))))

  keep <- mzs >= mz_range[1] & mzs <= mz_range[2]
  mzs <- mzs[keep]; mean_int <- mean_int[keep]

  # --- (1) mean-spectrum intensity of each m/z in the ion list
  res$mean_intensity <- vapply(res$mz, function(x) {
    w <- abs(mzs - x) <= x * ppm_tol * 1e-6
    if (any(w)) max(mean_int[w]) else NA_real_
  }, numeric(1))

  # --- (2) peak detection on the mean spectrum (local maxima > noise * S/N)
  noise <- median(abs(mean_int - median(mean_int))) * 1.4826
  is_pk <- c(FALSE, diff(sign(diff(mean_int))) == -2, FALSE) & mean_int > median(mean_int) + snr_min * noise
  pk <- data.frame(mz = mzs[is_pk], intensity = mean_int[is_pk])
  pk <- cbind(pk, mdf_classify(pk$mz)[, -1])
  pk$possible_DHB_cluster <- match_dhb(pk$mz)

  share <- function(df) {
    t <- tapply(df$intensity, df$mdf_class, sum)
    data.frame(Category = names(t), summed_intensity = as.numeric(t),
               percent_intensity = round(100 * as.numeric(t) / sum(t), 1),
               n_peaks = as.numeric(table(df$mdf_class)[names(t)]))
  }
  share_all <- share(pk)
  cat("\n==== imzML: intensity share (detected peaks, mean spectrum) ====\n")
  print(share_all, row.names = FALSE)

  xl <- res[!is.na(res$mean_intensity), ]
  xl$intensity <- xl$mean_intensity
  share_listed <- share(xl)
  cat("\n==== imzML: intensity share (ions of the input list) ====\n")
  print(share_listed, row.names = FALSE)

  p3 <- ggplot(pk, aes(mz, intensity, colour = mdf_class)) +
    geom_segment(aes(xend = mz, yend = 0)) +
    scale_colour_manual(values = cols, name = NULL) +
    labs(x = expression(italic(m/z)), y = "Mean intensity",
         title = "Mean spectrum - detected peaks coloured by mass defect class") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "mean_spectrum_mdf.png"), p3, width = 9, height = 4.5, dpi = 300)

  write_xlsx(list(Peaks_all = pk, Intensity_share_peaks = share_all,
                  Excel_list_intensity = xl, Intensity_share_listed = share_listed),
             file.path(out_dir, "mass_defect_imzML.xlsx"))
} else if (!is.na(imzml_file)) {
  warning("imzML file not found: ", imzml_file)
}

cat("\nDone. Results in:", normalizePath(out_dir), "\n")
