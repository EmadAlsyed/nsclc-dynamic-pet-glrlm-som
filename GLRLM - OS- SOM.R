###############################################################################
##  Dynamic 18F-FDG PET GLRLM radiomics -> Kohonen SOM -> overall survival
##  REVISED ANALYSIS (addresses Reviewers 1-3)
##
##  Inputs (same folder as this script):
##    Onlinedata.csv            survival (months) + 16 GLRLM features, 35 rows
##    pone.0239438.s002.xlsx    S1 Data of Noortman et al. 2020 (sheet "Clinical data")
##                              -> vital status (event) and, if present, patient ID
##
##  Everything the revised manuscript reports is written to NSCLC_SOM_Outputs/
##  (figures as 600-dpi PNG, tables as CSV, and manuscript_numbers.txt).
##
##  Section map (reviewer comment -> section)
##   1  Import, event status, patient vs tumour ........ R1.1, R2.1, R3
##   2  Main SOM (identical to submitted analysis) ...... reproducibility, R1.4, R1.5
##   3  SOM quality: QE, TE, weights vs data ............ R3
##   4  Stability: seeds, bootstrap, grid sizes ......... R2.2, R3
##   5  Choice of k and U-matrix ridge ................... R1.5, R2.3, R3
##   6  Benchmarks: PCA, k-means, hclust ................. R1 (overall), R3
##   7  Survival: KM, HR, PH test, RMST, permutation ..... R1.1, R1.6, R2.5, R2.6, R3
##   8  Node-level summaries (censoring-aware) ........... R1.7, R2.4
##   9  Figures (node numbers, grouped features) ......... R2 minor, R1.2
##  10  manuscript_numbers.txt + sessionInfo
###############################################################################

options(stringsAsFactors = FALSE, width = 120)

## ---- 0. Settings ------------------------------------------------------------
DATA_CSV   <- "Onlinedata.csv"
S1_XLSX    <- "pone.0239438.s002.xlsx"
S1_SHEET   <- "Clinical data"
OUT        <- "NSCLC_SOM_Outputs"
SEED       <- 2024            # seed of the submitted analysis
GX <- 3; GY <- 3; K <- 2      # grid and number of super-clusters (as submitted)
RLEN  <- 1000                 # training epochs (passes through the data)
ALPHA <- c(0.05, 0.01)        # learning rate, start -> end
N_SEEDS     <- 200            # random initialisations for the stability analysis
N_BOOT      <- 200            # bootstrap resamples for the stability analysis
N_SEED_GRID <- 30             # seeds per grid size
N_PERM      <- 10000          # permutations
N_BOOT_RMST <- 2000           # bootstrap resamples for the RMST difference
GRIDS <- list(c(2,2), c(2,3), c(3,3), c(3,4), c(4,4))

pkgs <- c("kohonen", "ggplot2", "RColorBrewer", "survival", "survminer",
          "readxl", "cluster")
new <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
dir.create(OUT, showWarnings = FALSE)
op <- function(f) file.path(OUT, f)
NUM <- character(0)                       # lines for manuscript_numbers.txt
note <- function(...) { s <- sprintf(...); NUM <<- c(NUM, s); cat(s, "\n") }
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
qs <- function(x) sprintf("median %.2f (IQR %.2f-%.2f; range %.2f-%.2f)",
                          median(x), quantile(x, .25), quantile(x, .75), min(x), max(x))

## ---- 1. Import data, vital status, patient identifier ------------------------
dat <- read.csv(DATA_CSV, check.names = FALSE)
names(dat) <- trimws(names(dat))
dat[] <- lapply(dat, function(x) as.numeric(as.character(x)))
dat <- dat[is.finite(dat$survival), ]
feat <- setdiff(names(dat), "survival")
stopifnot(length(feat) == 16, nrow(dat) == 35)

## Feature families (used to group features in the figures; reviewer 2, Fig. 3)
fam <- c(SRE = "Run emphasis", LRE = "Run emphasis",
         GLN = "Non-uniformity", GLNN = "Non-uniformity",
         RLN = "Non-uniformity", RLNN = "Non-uniformity",
         RP = "Run percentage",
         GLV = "Variance / entropy", RV = "Variance / entropy", RE = "Variance / entropy",
         LGLRE = "Gray-level emphasis", HGLRE = "Gray-level emphasis",
         SRLGLE = "Joint emphasis", SRHGLE = "Joint emphasis",
         LRLGLE = "Joint emphasis", LRHGLE = "Joint emphasis")
fam_levels <- c("Run emphasis", "Non-uniformity", "Run percentage",
                "Variance / entropy", "Gray-level emphasis", "Joint emphasis")
feat_ord <- feat[order(match(fam[feat], fam_levels))]
feat_names <- c(SRE = "Short run emphasis", LRE = "Long run emphasis",
                GLN = "Gray level non-uniformity", GLNN = "Gray level non-uniformity normalized",
                RLN = "Run length non-uniformity", RLNN = "Run length non-uniformity normalized",
                RP = "Run percentage", GLV = "Gray level variance", RV = "Run variance",
                RE = "Run entropy", LGLRE = "Low gray level run emphasis",
                HGLRE = "High gray level run emphasis",
                SRLGLE = "Short run low gray level emphasis",
                SRHGLE = "Short run high gray level emphasis",
                LRLGLE = "Long run low gray level emphasis",
                LRHGLE = "Long run high gray level emphasis")
write.csv(data.frame(Abbreviation = names(feat_names), Name = unname(feat_names),
                     Family = unname(fam[names(feat_names)])),
          op("Table_feature_abbreviations.csv"), row.names = FALSE)

## Vital status from S1 Data (row order; identical checks to the original script)
s1 <- read_excel(S1_XLSX, sheet = S1_SHEET)
s1 <- s1[!is.na(s1$`Overall survival`), ]
s1 <- s1[!is.na(suppressWarnings(as.integer(s1$`Overall survival`))), ]   # drop legend row
stopifnot(nrow(s1) == nrow(dat))
stopifnot(all(abs(as.numeric(s1$`Time overall survival (months)`) - dat$survival) < 1e-6))
dat$event <- as.integer(s1$`Overall survival`)

## Patient identifier (reviewer 2.1): OS is a patient-level endpoint.
## In S1 Data the first column of "Clinical data" holds the lesion ID; lesions of the same patient
## carry a letter suffix (26A, 26B). Stripping the suffix gives the patient ID (35 lesions, 34 patients).
lesion_id <- trimws(as.character(s1[[1]]))
stopifnot(length(lesion_id) == nrow(dat), !anyNA(lesion_id))
dat$patient <- sub("[A-Za-z]+$", "", lesion_id)
cat("Lesions sharing a patient ID:", paste(lesion_id[dat$patient %in% dat$patient[duplicated(dat$patient)]], collapse = ", "), "\n")
note("Tumours: %d | patients: %d | patients with >1 tumour: %d",
     nrow(dat), length(unique(dat$patient)), sum(table(dat$patient) > 1))
note("Events: %d (%.1f%%) | censored: %d (%.1f%%)", sum(dat$event), 100 * mean(dat$event),
     sum(dat$event == 0), 100 * mean(dat$event == 0))
note("Tumours with survival time = 0 months: %d (rows %s; events: %s)",
     sum(dat$survival == 0), paste(which(dat$survival == 0), collapse = ","),
     paste(dat$event[dat$survival == 0], collapse = ","))
write.csv(cbind(tumour = seq_len(nrow(dat)), dat[, c("patient", "survival", "event")]),
          op("Table_cohort_status.csv"), row.names = FALSE)

X <- as.matrix(dat[, feat])
zs <- function(M, center, scale) sweep(sweep(M, 2, center, "-"), 2, scale, "/")
Z  <- zs(X, colMeans(X), apply(X, 2, sd))            # identical to scale(X)

## ---- 2. Helper functions ----------------------------------------------------
train_som <- function(Xm, gx, gy, seed, rlen = RLEN, alpha = ALPHA) {
  Zm <- zs(Xm, colMeans(Xm), apply(Xm, 2, sd))
  set.seed(seed)
  som(Zm, grid = somgrid(xdim = gx, ydim = gy, topo = "hexagonal"),
      rlen = rlen, alpha = alpha, keep.data = TRUE)
}
get_cl <- function(m, k) {                           # complete-linkage on the codebook
  hc <- cutree(hclust(dist(m$codes[[1]]), method = "complete"), k = k)
  list(node = hc, tumour = hc[m$unit.classif])
}
ari <- function(a, b) {                              # adjusted Rand index
  tab <- table(a, b); n <- sum(tab); sc <- function(x) sum(choose(x, 2))
  sa <- sc(rowSums(tab)); sb <- sc(colSums(tab)); sab <- sc(tab)
  e <- sa * sb / choose(n, 2); den <- 0.5 * (sa + sb) - e
  if (abs(den) < 1e-12) 1 else (sab - e) / den
}
## quantization error (mean Euclidean distance to the BMU, in z-score units)
## and topographic error (share of tumours whose 1st and 2nd BMU are not adjacent)
qe_te <- function(m, Zm) {
  cw <- m$codes[[1]]; pts <- m$grid$pts
  r <- t(apply(Zm, 1, function(x) {
    d <- colSums((t(cw) - x)^2); o <- order(d)[1:2]
    c(sqrt(d[o[1]]), sqrt(sum((pts[o[1], ] - pts[o[2], ])^2)) > 1.05)
  }))
  c(QE = mean(r[, 1]), TE = mean(r[, 2]))
}
adj_pairs <- function(m) {                           # adjacent node pairs on the hex grid
  p <- m$grid$pts; d <- as.matrix(dist(p)); which(d > 0 & d < 1.05 & upper.tri(d), arr.ind = TRUE)
}

## ---- 3. Main SOM: identical to the submitted analysis ------------------------
som_main <- train_som(X, GX, GY, SEED)
cl_main  <- get_cl(som_main, K)
dat$Node <- som_main$unit.classif
dat$Cluster <- factor(cl_main$tumour)
counts <- as.integer(table(factor(dat$Node, levels = 1:(GX * GY))))
cat("Node counts:", counts, "\n")
if (!identical(counts, c(2L, 6L, 6L, 2L, 2L, 9L, 2L, 4L, 2L)) ||
    !identical(as.integer(cl_main$node), c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)))
  warning("The map differs from the one in the submitted manuscript (R/kohonen version?). ",
          "Report THESE results and update Tables 1-3 and Figures 3-8 accordingly.")
cat("kohonen neighbourhood function used:", som_main$neighbourhood.fct, "\n")   # report exactly this in Methods
note("kohonen %s | neighbourhood function: %s | rlen (epochs): %d | alpha %.2f->%.2f | R %s",
     as.character(packageVersion("kohonen")), as.character(som_main$neighbourhood.fct),
     RLEN, ALPHA[1], ALPHA[2], R.version.string)
ref_tumour <- cl_main$tumour
ref_node   <- dat$Node

qt <- qe_te(som_main, Z)
note("Main SOM: quantization error = %.3f z-units | topographic error = %.3f", qt["QE"], qt["TE"])
note("Free parameters: %d codebook weights vs %d data values (%d tumours x %d features); ratio %.2f",
     GX * GY * 16, 35 * 16, 35, 16, GX * GY * 16 / (35 * 16))

## Collinearity of the 16 features
R <- cor(X, method = "pearson"); diag(R) <- NA
note("Feature collinearity: %d of %d pairs have |r| >= 0.90; median |r| = %.2f",
     sum(abs(R[upper.tri(R)]) >= 0.9, na.rm = TRUE), 120, median(abs(R[upper.tri(R)]), na.rm = TRUE))
pca <- prcomp(Z)
ve  <- pca$sdev^2 / sum(pca$sdev^2)
n80 <- which(cumsum(ve) >= 0.80)[1]; n90 <- which(cumsum(ve) >= 0.90)[1]
note("PCA: PC1 %.1f%%, PC2 %.1f%%; %d PCs explain >=80%% and %d PCs explain >=90%% of variance",
     100 * ve[1], 100 * ve[2], n80, n90)

## Node z-score tables (Tables 1-2) and U-matrix values
node_z <- as.data.frame(som_main$codes[[1]]); node_z$Node <- 1:(GX * GY)
node_z$N <- counts; node_z$SuperCluster <- cl_main$node
write.csv(round(node_z[, c("Node", "N", "SuperCluster", feat_ord)], 3),
          op("Table1_2_node_zscores.csv"), row.names = FALSE)
ap <- adj_pairs(som_main); cw <- som_main$codes[[1]]
umat <- sapply(1:(GX * GY), function(i) {
  nb <- c(ap[ap[, 1] == i, 2], ap[ap[, 2] == i, 1]); mean(sqrt(rowSums((cw[nb, , drop = FALSE] - rep(cw[i, ], each = length(nb)))^2))) })
pair_d  <- sqrt(rowSums((cw[ap[, 1], ] - cw[ap[, 2], ])^2))
between <- cl_main$node[ap[, 1]] != cl_main$node[ap[, 2]]
note("U-matrix ridge: mean distance between adjacent nodes in DIFFERENT super-clusters = %.2f vs %.2f within super-clusters (ratio %.2f)",
     mean(pair_d[between]), mean(pair_d[!between]), mean(pair_d[between]) / mean(pair_d[!between]))

## ---- 4. Stability (reviewer 2.2, reviewer 3) --------------------------------
## 4a. Random initialisations, same grid
st_seed <- data.frame(); CL_seed <- matrix(NA_integer_, nrow(dat), N_SEEDS); ND_seed <- CL_seed
for (s in seq_len(N_SEEDS)) {
  m <- train_som(X, GX, GY, SEED + s); cl <- get_cl(m, K)
  q <- qe_te(m, Z)
  st_seed <- rbind(st_seed, data.frame(seed = SEED + s, ARI_cluster = ari(ref_tumour, cl$tumour),
                                       ARI_node = ari(ref_node, m$unit.classif), QE = unname(q["QE"]), TE = unname(q["TE"])))
  CL_seed[, s] <- cl$tumour; ND_seed[, s] <- m$unit.classif
}
write.csv(st_seed, op("TableS_stability_seeds.csv"), row.names = FALSE)
note("Seeds (n=%d): ARI of super-cluster partition vs reference: %s; share of seeds with ARI >= 0.80: %.0f%%; with ARI = 1: %.0f%%",
     N_SEEDS, qs(st_seed$ARI_cluster), 100 * mean(st_seed$ARI_cluster >= 0.8), 100 * mean(st_seed$ARI_cluster > 0.999))
note("Seeds: ARI of node assignment vs reference: %s", qs(st_seed$ARI_node))
note("Seeds: quantization error %s; topographic error %s", qs(st_seed$QE), qs(st_seed$TE))
## consensus (co-clustering frequency of tumours into the same super-cluster)
cons <- Reduce(`+`, lapply(seq_len(N_SEEDS), function(s) outer(CL_seed[, s], CL_seed[, s], "=="))) / N_SEEDS
ordr <- order(ref_tumour, dat$Node)
note("Consensus: mean co-clustering frequency of tumour pairs that share a super-cluster in the reference map = %.2f; of pairs that do not = %.2f",
     mean(cons[outer(ref_tumour, ref_tumour, "==") & upper.tri(cons)]),
     mean(cons[outer(ref_tumour, ref_tumour, "!=") & upper.tri(cons)]))

## 4b. Bootstrap (re-standardise and retrain on resampled tumours; project all 35 tumours)
st_boot <- data.frame()
for (b in seq_len(N_BOOT)) {
  set.seed(5000 + b); idx <- sample(nrow(X), replace = TRUE)
  Xb <- X[idx, ]; sdb <- apply(Xb, 2, sd); if (any(sdb == 0)) next
  m <- train_som(Xb, GX, GY, 6000 + b)
  Zall <- zs(X, colMeans(Xb), sdb)
  pr <- predict(m, newdata = Zall)$unit.classif
  cl <- cutree(hclust(dist(m$codes[[1]])), K)[pr]
  st_boot <- rbind(st_boot, data.frame(boot = b, ARI_cluster = ari(ref_tumour, cl)))
}
write.csv(st_boot, op("TableS_stability_bootstrap.csv"), row.names = FALSE)
note("Bootstrap (n=%d): ARI of super-cluster partition vs reference: %s; share with ARI >= 0.80: %.0f%%",
     nrow(st_boot), qs(st_boot$ARI_cluster), 100 * mean(st_boot$ARI_cluster >= 0.8))

## 4c. Grid size
LR_p <- function(d, g) { f <- survdiff(Surv(survival, event) ~ g, data = d); 1 - pchisq(f$chisq, length(f$n) - 1) }
grid_res <- data.frame()
for (g in GRIDS) {
  for (s in seq_len(N_SEED_GRID)) {
    m <- train_som(X, g[1], g[2], SEED + 100 + s); cl <- get_cl(m, K); q <- qe_te(m, Z)
    pp <- if (length(unique(cl$tumour)) == 2) LR_p(dat, factor(cl$tumour)) else NA
    grid_res <- rbind(grid_res, data.frame(grid = paste0(g[1], "x", g[2]), nodes = g[1] * g[2],
                                           seed = SEED + 100 + s, ARI_vs_3x3_reference = ari(ref_tumour, cl$tumour),
                                           QE = unname(q["QE"]), TE = unname(q["TE"]), logrank_p = pp,
                                           min_node_n = min(table(m$unit.classif)), empty_nodes = g[1] * g[2] - length(unique(m$unit.classif))))
  }
}
write.csv(grid_res, op("TableS_grid_size_sensitivity.csv"), row.names = FALSE)
grid_sum <- do.call(rbind, lapply(split(grid_res, grid_res$grid), function(d) data.frame(
  grid = d$grid[1], nodes = d$nodes[1], median_ARI = median(d$ARI_vs_3x3_reference),
  median_QE = median(d$QE), median_TE = median(d$TE), median_empty_nodes = median(d$empty_nodes),
  logrank_p_min = min(d$logrank_p, na.rm = TRUE), logrank_p_max = max(d$logrank_p, na.rm = TRUE))))
grid_sum <- grid_sum[order(grid_sum$nodes), ]
write.csv(grid_sum, op("TableS_grid_size_summary.csv"), row.names = FALSE)
print(grid_sum)
note("Grid-size sensitivity: log-rank p ranged from %.3f to %.3f across all grids and seeds (min over grids).",
     min(grid_res$logrank_p, na.rm = TRUE), max(grid_res$logrank_p, na.rm = TRUE))

## ---- 5. Choice of k (reviewer 1.5, 2.3, 3) -----------------------------------
dZ <- dist(Z)
k_tab <- do.call(rbind, lapply(2:5, function(k) {
  cl <- get_cl(som_main, k)$tumour
  data.frame(k = k, mean_silhouette = mean(silhouette(cl, dZ)[, "sil_width"]),
             smallest_cluster_n = min(table(cl)), cluster_sizes = paste(table(cl), collapse = "/"))
}))
write.csv(k_tab, op("TableS_choice_of_k.csv"), row.names = FALSE)
print(k_tab)
note("Silhouette (SOM super-clusters, z-score space): k=2: %.3f; k=3: %.3f; k=4: %.3f; k=5: %.3f",
     k_tab$mean_silhouette[1], k_tab$mean_silhouette[2], k_tab$mean_silhouette[3], k_tab$mean_silhouette[4])

## ---- 6. Benchmarks: PCA + k-means, k-means, hierarchical (reviewer 3) --------
set.seed(SEED)
bench <- list(
  "SOM + complete linkage (3x3, k=2)" = as.integer(ref_tumour),
  "k-means on z-scores (k=2)"          = kmeans(Z, 2, nstart = 100)$cluster,
  "PCA (>=80% var) + k-means (k=2)"    = kmeans(pca$x[, 1:n80, drop = FALSE], 2, nstart = 100)$cluster,
  "Hierarchical complete linkage on z-scores (k=2)" = cutree(hclust(dZ, "complete"), 2),
  "Hierarchical Ward.D2 on z-scores (k=2)"         = cutree(hclust(dZ, "ward.D2"), 2))

## ---- 7. Survival analysis (reviewers 1.1, 1.6, 2.1, 2.5, 2.6, 3) -------------
rmst1 <- function(t, e, tau) {                       # restricted mean survival time from the KM curve
  f <- survfit(Surv(t, e) ~ 1); i <- f$time < tau
  sum(diff(c(0, f$time[i], tau)) * c(1, f$surv[i]))
}
rmst_diff <- function(t, e, g, tau) {
  lv <- levels(g); rmst1(t[g == lv[2]], e[g == lv[2]], tau) - rmst1(t[g == lv[1]], e[g == lv[1]], tau)
}
perm_p <- function(d, stat, N = N_PERM) {            # permute (time, event) pairs over group labels
  obs <- stat(d$survival, d$event, d$g); n <- nrow(d); cnt <- 0
  for (i in seq_len(N)) {
    j <- sample.int(n); if (abs(stat(d$survival[j], d$event[j], d$g)) >= abs(obs) - 1e-12) cnt <- cnt + 1
  }
  (cnt + 1) / (N + 1)
}
one_analysis <- function(d, label, do_perm = TRUE, do_boot = TRUE) {
  d$g <- factor(d$g); stopifnot(nlevels(d$g) == 2)
  f  <- Surv(d$survival, d$event)
  lr <- survdiff(f ~ d$g); p_lr <- 1 - pchisq(lr$chisq, 1)
  pp <- survdiff(f ~ d$g, rho = 1); p_pp <- 1 - pchisq(pp$chisq, 1)
  cx <- coxph(f ~ d$g, ties = "efron"); ci <- summary(cx)$conf.int[1, ]
  zt <- cox.zph(cx)$table; p_ph <- zt[nrow(zt), "p"]
  cxr <- coxph(f ~ d$g + cluster(d$patient)); cir <- summary(cxr)$conf.int[1, ]
  tau <- min(tapply(d$survival, d$g, max))
  r_d <- rmst_diff(d$survival, d$event, d$g, tau)
  ci_r <- c(NA, NA)
  if (do_boot) {
    set.seed(SEED + 7)
    bs <- replicate(N_BOOT_RMST, {
      ib <- unlist(lapply(levels(d$g), function(l) sample(which(d$g == l), replace = TRUE)))
      rmst_diff(d$survival[ib], d$event[ib], d$g[ib], tau) })
    ci_r <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  }
  p_perm_lr <- p_perm_rmst <- NA
  if (do_perm) {
    set.seed(SEED + 11); p_perm_lr   <- perm_p(d, function(t, e, g) survdiff(Surv(t, e) ~ g)$chisq)
    set.seed(SEED + 13); p_perm_rmst <- perm_p(d, function(t, e, g) rmst_diff(t, e, g, tau))
  }
  data.frame(analysis = label, n = nrow(d), patients = length(unique(d$patient)), events = sum(d$event),
             events_g1 = sum(d$event[d$g == levels(d$g)[1]]), n_g1 = sum(d$g == levels(d$g)[1]),
             events_g2 = sum(d$event[d$g == levels(d$g)[2]]), n_g2 = sum(d$g == levels(d$g)[2]),
             logrank_chisq = lr$chisq, logrank_p = p_lr, petopeto_p = p_pp,
             HR_g2_vs_g1 = ci[1], HR_lo = ci[3], HR_hi = ci[4], cox_p = summary(cx)$coefficients[1, 5],
             HR_robust_lo = cir[3], HR_robust_hi = cir[4], schoenfeld_global_p = p_ph,
             tau_months = tau, RMST_diff_months = r_d, RMST_lo = unname(ci_r[1]), RMST_hi = unname(ci_r[2]),
             perm_p_logrank = p_perm_lr, perm_p_RMST = p_perm_rmst)
}

## One tumour per patient: every possible choice when a patient has >1 tumour
one_per_patient_sets <- function(d) {
  ids <- unique(d$patient[duplicated(d$patient)]); base <- which(!d$patient %in% ids)
  if (!length(ids)) return(list(seq_len(nrow(d))))
  rows <- lapply(ids, function(i) which(d$patient == i)); cmb <- expand.grid(lapply(rows, seq_along))
  lapply(seq_len(nrow(cmb)), function(r)
    sort(c(base, mapply(function(rr, j) rr[j], rows, as.integer(unlist(cmb[r, ]))))))
}
dat$g <- dat$Cluster
sets <- one_per_patient_sets(dat)
res <- list()
res[["A. all tumours (patient-clustered robust SE in Cox)"]] <- one_analysis(dat, "A. all tumours (35)")
for (i in seq_along(sets))
  res[[paste0("B", i)]] <- one_analysis(dat[sets[[i]], ], sprintf("B%d. one tumour per patient (choice %d)", i, i))
res[["C. one tumour per patient, excluding survival = 0"]] <-
  one_analysis(dat[sets[[1]], ][dat$survival[sets[[1]]] > 0, ], "C. one tumour/patient, survival > 0")
surv_tab <- do.call(rbind, res); rownames(surv_tab) <- NULL
write.csv(surv_tab, op("Table4_survival_analyses.csv"), row.names = FALSE)
print(t(surv_tab))
prim <- surv_tab[2, ]                                   # primary: one tumour per patient (choice 1)
if (nrow(surv_tab) == 3) prim <- surv_tab[1, ]          # no duplicate patient found: all tumours = patients
note("PRIMARY (%s): n=%d, events %d; log-rank chi2=%.2f, p=%s; HR(SC2 vs SC1)=%.2f (95%% CI %.2f-%.2f); Schoenfeld global p=%s",
     prim$analysis, prim$n, prim$events, prim$logrank_chisq, fmtp(prim$logrank_p),
     prim$HR_g2_vs_g1, prim$HR_lo, prim$HR_hi, fmtp(prim$schoenfeld_global_p))
note("PRIMARY: RMST difference (SC2 - SC1) up to tau=%.0f months = %.1f months (95%% bootstrap CI %.1f to %.1f); permutation p (log-rank) = %s; permutation p (RMST) = %s; Peto-Peto p = %s",
     prim$tau_months, prim$RMST_diff_months, prim$RMST_lo, prim$RMST_hi,
     fmtp(prim$perm_p_logrank), fmtp(prim$perm_p_RMST), fmtp(prim$petopeto_p))
note("All %d one-tumour-per-patient choices: log-rank p range %s-%s; HR range %.2f-%.2f",
     sum(grepl("^B", surv_tab$analysis)),
     fmtp(min(surv_tab$logrank_p[grepl("^B", surv_tab$analysis)])), fmtp(max(surv_tab$logrank_p[grepl("^B", surv_tab$analysis)])),
     min(surv_tab$HR_g2_vs_g1[grepl("^B", surv_tab$analysis)]), max(surv_tab$HR_g2_vs_g1[grepl("^B", surv_tab$analysis)]))
## sensitivity at tau = 60 months
d_prim <- dat[sets[[1]], ]; d_prim$g <- factor(d_prim$g)
set.seed(SEED + 17)
bs60 <- replicate(N_BOOT_RMST, { ib <- unlist(lapply(levels(d_prim$g), function(l) sample(which(d_prim$g == l), replace = TRUE)))
rmst_diff(d_prim$survival[ib], d_prim$event[ib], d_prim$g[ib], 60) })
note("RMST difference at tau = 60 months = %.1f (95%% CI %.1f to %.1f)",
     rmst_diff(d_prim$survival, d_prim$event, d_prim$g, 60), quantile(bs60, .025), quantile(bs60, .975))

## Kaplan-Meier tables (primary patient-level set, one lesion per patient): medians with CI, censoring, risk table
dkm <- dat[sets[[1]], ]
fit_all <- survfit(Surv(survival, event) ~ Cluster, data = dkm)
km_tab <- summary(fit_all)$table
write.csv(as.data.frame(km_tab), op("Table_KM_summary_by_supercluster.csv"))
print(km_tab)
steps <- summary(fit_all, censored = TRUE)
write.csv(data.frame(stratum = steps$strata, time = steps$time, n.risk = steps$n.risk,
                     n.event = steps$n.event, n.censor = steps$n.censor, surv = steps$surv,
                     lower = steps$lower, upper = steps$upper),
          op("Table_KM_steps_by_supercluster.csv"), row.names = FALSE)
for (l in levels(dkm$Cluster)) {
  dd <- dkm[dkm$Cluster == l, ]
  f <- survfit(Surv(survival, event) ~ 1, data = dd); tb <- summary(f)$table
  rf <- survfit(Surv(survival, 1 - event) ~ 1, data = dd)               # reverse KM: median follow-up
  note("Primary set, super-cluster %s: n=%d, events=%d, censored=%d (%.0f%%); KM median = %s months (95%% CI %s-%s); lowest KM estimate = %.3f; median follow-up (reverse KM) = %s months",
       l, nrow(dd), sum(dd$event), sum(dd$event == 0), 100 * mean(dd$event == 0),
       ifelse(is.na(tb["median"]), "not reached", sprintf("%.0f", tb["median"])),
       ifelse(is.na(tb["0.95LCL"]), "NA", sprintf("%.0f", tb["0.95LCL"])),
       ifelse(is.na(tb["0.95UCL"]), "NA", sprintf("%.0f", tb["0.95UCL"])),
       min(f$surv), ifelse(is.na(summary(rf)$table["median"]), "NA", sprintf("%.0f", summary(rf)$table["median"])))
}
## when do the curves cross?
sf <- lapply(levels(dkm$Cluster), function(l) { f <- survfit(Surv(survival, event) ~ 1, data = dkm[dkm$Cluster == l, ]); stepfun(f$time, c(1, f$surv)) })
tg <- seq(0, max(dkm$survival), by = 0.5); dd <- sf[[2]](tg) - sf[[1]](tg)
nz <- which(dd != 0); sg <- sign(dd[nz])
cr <- tg[nz][which(diff(sg) != 0) + 1]               # first grid time after each sign change
note("Kaplan-Meier curves cross at approximately %s months", ifelse(length(cr), paste(round(cr, 1), collapse = ", "), "(no crossing)"))

## ---- 6 (cont.) benchmark table with survival results -------------------------
bt <- do.call(rbind, lapply(names(bench), function(nm) {
  cl <- bench[[nm]]; d <- dat; d$g <- factor(cl)
  cx <- coxph(Surv(survival, event) ~ g, data = d); ci <- summary(cx)$conf.int[1, ]
  data.frame(method = nm, sizes = paste(table(cl), collapse = "/"), ARI_vs_SOM = ari(ref_tumour, cl),
             mean_silhouette = mean(silhouette(as.integer(cl), dZ)[, "sil_width"]),
             logrank_p = LR_p(d, d$g), HR = ci[1], HR_lo = ci[3], HR_hi = ci[4])
}))
write.csv(bt, op("Table5_benchmark_methods.csv"), row.names = FALSE)
print(bt)
note("Benchmarks: ARI vs SOM partition: k-means %.2f; PCA+k-means %.2f; hierarchical complete %.2f; Ward.D2 %.2f. Log-rank p: %s",
     bt$ARI_vs_SOM[2], bt$ARI_vs_SOM[3], bt$ARI_vs_SOM[4], bt$ARI_vs_SOM[5], paste(fmtp(bt$logrank_p), collapse = " / "))

## ---- 8. Node-level summaries (censoring-aware; reviewers 1.7, 2.4) -----------
node_tab <- do.call(rbind, lapply(1:(GX * GY), function(i) {
  d <- dat[dat$Node == i, ]
  data.frame(Node = i, SuperCluster = cl_main$node[i], N = nrow(d), Events = sum(d$event),
             Censored = sum(d$event == 0), Pct_censored = round(100 * mean(d$event == 0)),
             Mean_observed_time = round(mean(d$survival), 1),        # event OR censoring time; NOT survival
             Median_observed_time = median(d$survival))
}))
write.csv(node_tab, op("Table3_node_summary_censoring_aware.csv"), row.names = FALSE)
print(node_tab)

## ---- 9. Figures --------------------------------------------------------------
pts <- som_main$grid$pts
hex <- function(vals) {
  r <- 1 / sqrt(3); a <- pi / 6 + (0:5) * pi / 3
  do.call(rbind, lapply(seq_len(nrow(pts)), function(i)
    data.frame(node = i, x = pts[i, 1] + r * cos(a), y = pts[i, 2] + r * sin(a), val = vals[i])))
}
ctr <- data.frame(node = 1:nrow(pts), x = pts[, 1], y = pts[, 2])
save_gg <- function(p, f, w, h) ggsave(op(f), p, width = w, height = h, dpi = 600, bg = "white")
th <- theme_void(base_size = 11) + theme(plot.title = element_text(face = "bold", hjust = .5))
node_lab <- sprintf("N%d\nn=%d", 1:(GX * GY), counts)

## Fig. 2  training progress (convergence only; not a stability measure)
ch <- data.frame(epoch = seq_along(som_main$changes[, 1]), dist = som_main$changes[, 1])
save_gg(ggplot(ch, aes(epoch, dist)) + geom_line(linewidth = .6) + theme_bw(base_size = 11) +
          labs(x = "Training epoch", y = "Mean distance to BMU (kohonen units)",
               title = "SOM training progress (convergence)"), "Fig2_SOM_training.png", 6, 4.2)

## Fig. 3  codebook heat map, features grouped by family (replaces the 16-colour pie plot)
cb <- do.call(rbind, lapply(feat_ord, function(f) data.frame(
  node = factor(paste0("N", 1:(GX * GY)), levels = paste0("N", 1:(GX * GY))), feature = factor(f, levels = feat_ord),
  family = factor(fam[f], levels = fam_levels), z = som_main$codes[[1]][, f])))
save_gg(ggplot(cb, aes(feature, node, fill = z)) + geom_tile(colour = "white") +
          geom_text(aes(label = sprintf("%.1f", z)), size = 2.2) +
          scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-3.7, 3.7), name = "z-score") +
          facet_grid(. ~ family, scales = "free_x", space = "free_x", labeller = label_wrap_gen(10)) + scale_y_discrete(limits = rev) +
          labs(x = NULL, y = NULL, title = "Codebook vectors (z-scores) by node, features grouped by family") +
          theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1),
                                                strip.text = element_text(size = 7), panel.grid = element_blank()),
        "Fig3_codebook_heatmap.png", 11, 5)

## Fig. 4  tumours per node
save_gg(ggplot() + geom_polygon(data = hex(counts), aes(x, y, group = node, fill = val), colour = "grey30") +
          geom_text(data = cbind(ctr, lab = node_lab), aes(x, y, label = lab), size = 3) +
          scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Tumors") +
          coord_equal() + labs(title = "Tumors per SOM node") + th, "Fig4_counts.png", 6, 5.2)

## Fig. 5  component planes (one panel per feature, shared diverging scale)
cp <- do.call(rbind, lapply(feat_ord, function(f) { h <- hex(som_main$codes[[1]][, f]); h$feature <- factor(f, levels = feat_ord); h }))
save_gg(ggplot(cp, aes(x, y, group = interaction(node, feature), fill = val)) + geom_polygon(colour = "grey40", linewidth = .2) +
          scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "z-score") +
          facet_wrap(~feature, ncol = 4) + coord_equal() +
          labs(title = "Component planes (node codebook z-scores)") + th, "Fig5_component_planes.png", 9, 9)

## Fig. 6  U-matrix
save_gg(ggplot() + geom_polygon(data = hex(umat), aes(x, y, group = node, fill = val), colour = "grey30") +
          geom_text(data = cbind(ctr, lab = sprintf("N%d", 1:(GX * GY))), aes(x, y, label = lab), size = 3) +
          scale_fill_distiller(palette = "Greys", direction = 1, name = "Mean distance\nto neighbours") +
          coord_equal() + labs(title = "U-matrix") + th, "Fig6_Umatrix.png", 6, 5.2)

## Fig. 7 (supplement)  observed time per node: NOT survival (20/35 censored)
save_gg(ggplot() + geom_polygon(data = hex(node_tab$Mean_observed_time), aes(x, y, group = node, fill = val), colour = "grey30") +
          geom_text(data = cbind(ctr, lab = sprintf("N%d\n%d/%d events", 1:(GX * GY), node_tab$Events, node_tab$N)), aes(x, y, label = lab), size = 3) +
          scale_fill_distiller(palette = "Blues", direction = 1, name = "Mean observed\ntime (months)") +
          coord_equal() + labs(title = "Mean observed time per node (descriptive)") + th,
        "FigS_observed_time_map.png", 6, 5.2)

## Fig. 8  super-clusters with node numbers
sc_col <- c("1" = "#66C2A5", "2" = "#FC8D62")
h8 <- hex(factor(cl_main$node));
save_gg(ggplot() + geom_polygon(data = h8, aes(x, y, group = node, fill = val), colour = "grey30") +
          geom_text(data = cbind(ctr, lab = node_lab), aes(x, y, label = lab), size = 3.2) +
          scale_fill_manual(values = sc_col, name = "Super-cluster") + coord_equal() +
          labs(title = sprintf("Super-clusters: SC1 = nodes %s; SC2 = nodes %s",
                               paste(which(cl_main$node == 1), collapse = ","), paste(which(cl_main$node == 2), collapse = ","))) + th,
        "Fig8_superclusters.png", 6.4, 5.2)

## Fig. 9  Kaplan-Meier (one tumour per patient shown as primary; all tumours in the supplement)
mk_km <- function(d, ptxt, f) {
  d$g <- factor(d$g, labels = c("Super-cluster 1", "Super-cluster 2"))
  fit <- survfit(Surv(survival, event) ~ g, data = d)
  p <- ggsurvplot(fit, data = d, pval = ptxt, pval.size = 3.4, risk.table = TRUE, conf.int = TRUE,
                  palette = unname(sc_col), censor = TRUE, surv.median.line = "hv",
                  xlab = "Overall survival (months)", legend.title = "", legend.labs = levels(d$g),
                  break.time.by = 25, ggtheme = theme_bw(), risk.table.height = .25)
  png(op(f), width = 7, height = 7.5, units = "in", res = 600); print(p); dev.off()
}
txt <- sprintf("Log-rank p = %s\nHR (SC2 vs SC1) = %.2f (95%% CI %.2f-%.2f)\nPermutation p = %s; Schoenfeld p = %s\nRMST diff. = %.1f mo (95%% CI %.1f to %.1f)",
               fmtp(prim$logrank_p), prim$HR_g2_vs_g1, prim$HR_lo, prim$HR_hi, fmtp(prim$perm_p_logrank),
               fmtp(prim$schoenfeld_global_p), prim$RMST_diff_months, prim$RMST_lo, prim$RMST_hi)
mk_km(dat[sets[[1]], ], txt, "Fig9_KM_primary.png")
mk_km(dat, sprintf("Log-rank p = %s (all %d lesions)", fmtp(surv_tab$logrank_p[1]), nrow(dat)), "FigS_KM_all_tumours.png")

## Supplementary figures
sup <- rbind(data.frame(analysis = "Random initializations", ARI = st_seed$ARI_cluster),
             data.frame(analysis = "Bootstrap resamples", ARI = st_boot$ARI_cluster))
save_gg(ggplot(sup, aes(ARI)) + geom_histogram(bins = 25, fill = "grey60", colour = "white") + facet_wrap(~analysis, ncol = 1) +
          theme_bw(base_size = 11) + labs(x = "Adjusted Rand index vs reference partition", y = "Count",
                                          title = "Stability of the two-super-cluster partition"), "FigS_stability_ARI.png", 6, 6)
cm <- expand.grid(i = seq_along(ordr), j = seq_along(ordr)); cm$v <- as.vector(cons[ordr, ordr])
save_gg(ggplot(cm, aes(i, j, fill = v)) + geom_tile() + scale_fill_viridis_c(name = "Co-clustering\nfrequency") +
          coord_equal() + theme_minimal() + labs(x = "Lesions (ordered by reference super-cluster, node)", y = NULL,
                                                 title = "Consensus over random initializations"), "FigS_consensus.png", 6.5, 5.5)
gl <- grid_res; gl$grid <- factor(gl$grid, levels = unique(gl$grid[order(gl$nodes)]))
save_gg(ggplot(gl, aes(grid, ARI_vs_3x3_reference)) + geom_boxplot(fill = "grey85") + theme_bw(base_size = 11) +
          labs(x = "SOM grid", y = "ARI of 2-cluster partition vs reference", title = "Sensitivity to grid size"), "FigS_grid_size_ARI.png", 6, 4.2)
kp <- data.frame(k = k_tab$k, sil = k_tab$mean_silhouette)
save_gg(ggplot(kp, aes(k, sil)) + geom_line() + geom_point(size = 2) + theme_bw(base_size = 11) +
          labs(x = "Number of super-clusters k", y = "Mean silhouette width", title = "Choice of k"), "FigS_choice_of_k.png", 5, 4)
pc <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], SC = factor(ref_tumour), Node = factor(dat$Node))
save_gg(ggplot(pc, aes(PC1, PC2, colour = SC)) + geom_point(size = 2.2) + scale_colour_manual(values = sc_col, name = "SOM super-cluster") +
          theme_bw(base_size = 11) + labs(x = sprintf("PC1 (%.1f%%)", 100 * ve[1]), y = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
                                          title = "PCA of the 16 z-scored features"), "FigS_PCA.png", 6, 4.5)
Rm <- cor(X); rd <- expand.grid(a = factor(feat_ord, levels = feat_ord), b = factor(feat_ord, levels = feat_ord))
rd$r <- as.vector(Rm[feat_ord, feat_ord])
save_gg(ggplot(rd, aes(a, b, fill = r)) + geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1)) +
          theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + labs(x = NULL, y = NULL, fill = "Pearson r",
                                                                                                          title = "Collinearity of the 16 dynamic GLRLM features"), "FigS_correlation.png", 7, 6)

## ---- 10. Outputs ---------------------------------------------------------------
writeLines(NUM, op("manuscript_numbers.txt"))
writeLines(capture.output(sessionInfo()), op("sessionInfo.txt"))
cat("\nDone. Open", op("manuscript_numbers.txt"), "and copy the values into the revised manuscript.\n")