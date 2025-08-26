#Three‑level meta‑analysis of LDL thresholds (≥5 studies) – SE + Prediction I‑95
library(metafor)
setwd("C:/Users/ey270/Downloads/SoRS_meta_analysis_package")
## ---------------------------------------------------------------------------
## 1. Load data
## ---------------------------------------------------------------------------
df <- read.csv("StatinMeta-analysis.csv", check.names = TRUE)

names(df) <- make.names(names(df))      # valid column names
df$Goal   <- as.character(df$Goal)

## ---------------------------------------------------------------------------
## 2. Extract numeric LDL threshold from 'Goal'
## ---------------------------------------------------------------------------
extract_threshold <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  m <- regmatches(x, regexpr("[0-9]+\\.?[0-9]*", x))
  if (length(m) == 0) return(NA_real_)
  as.numeric(m)
}
df$Goal_threshold <- vapply(df$Goal, extract_threshold, numeric(1))

## ---------------------------------------------------------------------------
## 3. Build consistent study IDs
## ---------------------------------------------------------------------------
clean_study_id <- function(author) {
  a <- tolower(author)
  a <- gsub("[^a-z0-9]", "", a)     # remove punctuation / spaces
  a <- gsub("etal", "", a)          # drop “et al”
  a <- gsub("[0-9]{4}$", "", a)     # drop trailing year
  a
}
df$study_id <- vapply(df$author, clean_study_id, character(1))

## ---------------------------------------------------------------------------
## 4. Keep thresholds reported by ≥5 studies
## ---------------------------------------------------------------------------
df <- df[!is.na(df$Goal_threshold), ]
uniq_combo <- unique(data.frame(thr  = df$Goal_threshold,
                                stud = df$study_id,
                                stringsAsFactors = FALSE))
thr_counts    <- table(uniq_combo$thr)
eligible_thr  <- as.numeric(names(thr_counts[thr_counts >= 5]))

if (length(eligible_thr) == 0L) stop("No threshold is used in ≥5 studies.")
cat("Thresholds selected (≥5 studies):",
    paste(eligible_thr, collapse = ", "), "\n")

df <- df[df$Goal_threshold %in% eligible_thr, ]

## ---------------------------------------------------------------------------
## 5. Helper
## ---------------------------------------------------------------------------
inv_logit <- function(x) exp(x) / (1 + exp(x))

## ---------------------------------------------------------------------------
## 6. Meta‑analysis loop
## ---------------------------------------------------------------------------
heterogeneity_log <- data.frame()

for (thr in sort(unique(df$Goal_threshold))) {
  message("\n===== analysing threshold ", thr, " =====")
  tmp <- df[df$Goal_threshold == thr, ]
  # Track total participants for the current threshold
  total_participants <- sum(as.numeric(tmp$samplesize[ok]), na.rm = TRUE)
  cat(sprintf("  total participants for %s mg/dL threshold: %d\n", thr, total_participants))
  
  ## numeric vectors
  xi  <- suppressWarnings(as.numeric(tmp$numberofsors))
  ni  <- suppressWarnings(as.numeric(tmp$samplesize))
  lab <- paste(tmp$author)
  lab2<- paste(tmp$statins, tmp$Goal_details)          # forest‑plot labels
  
  ## validity filter
  ok  <- !is.na(xi) & !is.na(ni) & xi >= 0 & ni > 0 & xi <= ni
  xi  <- xi[ok]; ni <- ni[ok]; lab <- lab[ok]
  
  ## need ≥5 unique studies after filtering
  if (length(unique(tmp$study_id[ok])) < 5L) {
    message("  skipped – fewer than 5 studies with valid data")
    next
  }
  
  ## effect size calculation
  esc <- escalc(measure = "PLO", xi = xi, ni = ni)
  esc$study_id  <- factor(tmp$study_id[ok])
  esc$effect_id <- seq_along(xi)
  
  ## three‑ and two‑level REML models
  m3 <- rma.mv(yi, vi,
               random = list(~1 | effect_id, ~1 | study_id),
               data   = esc, method = "REML")
  m2 <- rma.mv(yi, vi,
               random = ~1 | study_id,
               data   = esc, method = "REML")
  
  ## model fit
  AIC_2 <- AIC(m2);  AIC_3 <- AIC(m3)
  BIC_2 <- BIC(m2);  BIC_3 <- BIC(m3)
  dAIC  <- AIC_2 - AIC_3
  dBIC  <- BIC_2 - BIC_3
  
  ## pooled prevalence, SE, CI, prediction interval
  prop     <- inv_logit(m3$b)                       # point estimate
  propCI   <- inv_logit(c(m3$ci.lb, m3$ci.ub))      # confidence interval
  se_logit <- sqrt(m3$vb[1, 1])                     # SE on logit scale
  se_prop  <- se_logit * prop * (1 - prop)          # delta‑method SE
  pred     <- predict(m3)                           # PI on logit scale
  pred_int <- inv_logit(c(pred$pi.lb, pred$pi.ub))  # PI on prevalence scale
  
  cat(sprintf("  pooled prevalence %.3f (SE %.4f) 95%% CI %.3f–%.3f\n",
              prop, se_prop, propCI[1], propCI[2]))
  cat(sprintf("  prediction interval %.3f–%.3f\n", pred_int[1], pred_int[2]))
  cat(sprintf("  AIC 2‑level %.2f | 3‑level %.2f (Δ %.2f)  |  ",
              AIC_2, AIC_3, dAIC))
  cat(sprintf("BIC 2‑level %.2f | 3‑level %.2f (Δ %.2f)\n",
              BIC_2, BIC_3, dBIC))
  
  ## heterogeneity
  tau2_w <- m3$sigma2[1]; tau2_b <- m3$sigma2[2]
  n_eff  <- nrow(esc)
  inv_v  <- 1/esc$vi
  typ_sv <- ((n_eff - 1) * sum(inv_v)) /
    (sum(inv_v)^2 - sum(inv_v^2))
  tot_v  <- typ_sv + tau2_w + tau2_b
  I2_w   <- 100 * tau2_w / tot_v
  I2_b   <- 100 * tau2_b / tot_v
  
  ## log table
  heterogeneity_log <- rbind(heterogeneity_log, data.frame(
    Threshold      = thr,
    Pooled         = prop,
    SE_pooled      = se_prop,
    CI_low         = propCI[1],
    CI_high        = propCI[2],
    PI_low         = pred_int[1],
    PI_high        = pred_int[2],
    Tau2_within    = tau2_w,
    Tau2_between   = tau2_b,
    I2_within      = I2_w,
    I2_between     = I2_b,
    AIC_2level     = AIC_2,
    AIC_3level     = AIC_3,
    dAIC           = dAIC,
    BIC_2level     = BIC_2,
    BIC_3level     = BIC_3,
    dBIC           = dBIC,
    stringsAsFactors = FALSE
  ))
  
  ## forest plot
  ticks_p <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.99)
  ticks_l <- log(ticks_p / (1 - ticks_p))
  pdf(paste0("Forest_", gsub("\\.", "_", thr), ".pdf"), 9, 7)
  
  par(mar = c(5, 4, 5, 2))
  metafor::forest(m3, slab = lab, ilab=lab2,
                  atransf = inv_logit,
                  at      = ticks_l,
                  xlab    = "Prevalence of sub‑optimal response",
                  xlim = c(-13, 5),
                  cex     = 0.8)
  axis(1, at = ticks_l,
       labels = formatC(ticks_p, format = "f", digits = 2), cex.axis = 0.8)
  title(main = paste("3‑level meta‑analysis – LDL goal", thr), cex.main = 1.1)
  
  dev.off()
  if (dev.cur() != 1) dev.off()
}

## ---------------------------------------------------------------------------
## 7. Save summary
## ---------------------------------------------------------------------------
write.csv(heterogeneity_log, "heterogeneity_summary_AIC_BIC.csv", row.names = FALSE)
cat("\n✓ heterogeneity_summary_AIC_BIC.csv written\n")
print(heterogeneity_log)

# Combined total across included thresholds
included_thresholds <- sort(unique(df$Goal_threshold))
total_combined_participants <- sum(as.numeric(df$samplesize[df$Goal_threshold %in% included_thresholds]), na.rm = TRUE)
cat(sprintf("\nTotal combined participants across thresholds %s: %d\n",
            paste(included_thresholds, collapse = " & "),
            total_combined_participants))










































#######Sensitivity analysis removing 0 values#################

###############################################################################
# Three‑level meta‑analysis of LDL thresholds (≥5 studies) – SE + Prediction I‑95
# Column names: samplesize, numberofsors, author, Goal
###############################################################################
library(metafor)

## ---------------------------------------------------------------------------
## 1. Load data
## ---------------------------------------------------------------------------
df <- read.csv("data_3_7.csv",
               stringsAsFactors = FALSE, check.names = TRUE)

names(df) <- make.names(names(df))      # valid column names
df$Goal   <- as.character(df$Goal)

## ---------------------------------------------------------------------------
## 2. Extract numeric LDL threshold from 'Goal'
## ---------------------------------------------------------------------------
extract_threshold <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  m <- regmatches(x, regexpr("[0-9]+\\.?[0-9]*", x))
  if (length(m) == 0) return(NA_real_)
  as.numeric(m)
}
df$Goal_threshold <- vapply(df$Goal, extract_threshold, numeric(1))

## ---------------------------------------------------------------------------
## 3. Build consistent study IDs
## ---------------------------------------------------------------------------
clean_study_id <- function(author) {
  a <- tolower(author)
  a <- gsub("[^a-z0-9]", "", a)     # remove punctuation / spaces
  a <- gsub("etal", "", a)          # drop “et al”
  a <- gsub("[0-9]{4}$", "", a)     # drop trailing year
  a
}
df$study_id <- vapply(df$author, clean_study_id, character(1))

## ---------------------------------------------------------------------------
## 4. Keep thresholds reported by ≥5 studies
## ---------------------------------------------------------------------------
df <- df[!is.na(df$Goal_threshold), ]
uniq_combo <- unique(data.frame(thr  = df$Goal_threshold,
                                stud = df$study_id,
                                stringsAsFactors = FALSE))
thr_counts    <- table(uniq_combo$thr)
eligible_thr  <- as.numeric(names(thr_counts[thr_counts >= 5]))

if (length(eligible_thr) == 0L) stop("No threshold is used in ≥5 studies.")
cat("Thresholds selected (≥5 studies):",
    paste(eligible_thr, collapse = ", "), "\n")

df <- df[df$Goal_threshold %in% eligible_thr, ]

## ---------------------------------------------------------------------------
## 5. Helper
## ---------------------------------------------------------------------------
inv_logit <- function(x) exp(x) / (1 + exp(x))

## --------------------------------------------------------------------------

























































heterogeneity_log_sens <- data.frame()  # To store sensitivity results across thresholds

for (thr in sort(unique(df$Goal_threshold))) {
  message("\n===== analysing threshold ", thr, " =====")
  tmp <- df[df$Goal_threshold == thr, ]
  
  ## Main meta-analysis (all data) -------------------
  xi  <- suppressWarnings(as.numeric(tmp$numberofsors))
  ni  <- suppressWarnings(as.numeric(tmp$samplesize))
  lab <- paste(tmp$author)
  
  ok  <- !is.na(xi) & !is.na(ni) & xi >= 0 & ni > 0 & xi <= ni
  xi  <- xi[ok]; ni <- ni[ok]; lab <- lab[ok]
  
  if (length(unique(tmp$study_id[ok])) < 5L) {
    message("  skipped – fewer than 5 studies with valid data")
    next
  }
  
  esc <- escalc(measure = "PLO", xi = xi, ni = ni)
  esc$study_id  <- factor(tmp$study_id[ok])
  esc$effect_id <- seq_along(xi)
  
  m3 <- rma.mv(yi, vi,
               random = list(~1 | effect_id, ~1 | study_id),
               data   = esc, method = "REML")
  m2 <- rma.mv(yi, vi,
               random = ~1 | study_id,
               data   = esc, method = "REML")
  
  AIC_2 <- AIC(m2);  AIC_3 <- AIC(m3)
  BIC_2 <- BIC(m2);  BIC_3 <- BIC(m3)
  dAIC  <- AIC_2 - AIC_3
  dBIC  <- BIC_2 - BIC_3
  
  prop     <- inv_logit(m3$b)
  propCI   <- inv_logit(c(m3$ci.lb, m3$ci.ub))
  se_logit <- sqrt(m3$vb[1, 1])
  se_prop  <- se_logit * prop * (1 - prop)
  pred     <- predict(m3)
  pred_int <- inv_logit(c(pred$pi.lb, pred$pi.ub))
  
  cat(sprintf("  pooled prevalence %.3f (SE %.4f) 95%% CI %.3f–%.3f\n",
              prop, se_prop, propCI[1], propCI[2]))
  cat(sprintf("  prediction interval %.3f–%.3f\n", pred_int[1], pred_int[2]))
  cat(sprintf("  AIC 2-level %.2f | 3-level %.2f (Δ %.2f)  |  ",
              AIC_2, AIC_3, dAIC))
  cat(sprintf("BIC 2-level %.2f | 3-level %.2f (Δ %.2f)\n",
              BIC_2, BIC_3, dBIC))
  
  tau2_w <- m3$sigma2[1]; tau2_b <- m3$sigma2[2]
  n_eff  <- nrow(esc)
  inv_v  <- 1/esc$vi
  typ_sv <- ((n_eff - 1) * sum(inv_v)) /
    (sum(inv_v)^2 - sum(inv_v^2))
  tot_v  <- typ_sv + tau2_w + tau2_b
  I2_w   <- 100 * tau2_w / tot_v
  I2_b   <- 100 * tau2_b / tot_v
  
  heterogeneity_log <- rbind(heterogeneity_log, data.frame(
    Threshold      = thr,
    Pooled         = prop,
    SE_pooled      = se_prop,
    CI_low         = propCI[1],
    CI_high        = propCI[2],
    PI_low         = pred_int[1],
    PI_high        = pred_int[2],
    Tau2_within    = tau2_w,
    Tau2_between   = tau2_b,
    I2_within      = I2_w,
    I2_between     = I2_b,
    AIC_2level     = AIC_2,
    AIC_3level     = AIC_3,
    dAIC           = dAIC,
    BIC_2level     = BIC_2,
    BIC_3level     = BIC_3,
    dBIC           = dBIC,
    stringsAsFactors = FALSE
  ))
  
  ## Sensitivity analysis: exclude numberofsors == 0 -----------
  tmp_sens <- tmp[tmp$numberofsors != 0, ]
  xi_sens  <- suppressWarnings(as.numeric(tmp_sens$numberofsors))
  ni_sens  <- suppressWarnings(as.numeric(tmp_sens$samplesize))
  lab_sens <- paste(tmp_sens$author)
  
  ok_sens <- !is.na(xi_sens) & !is.na(ni_sens) & xi_sens > 0 & ni_sens > 0 & xi_sens <= ni_sens
  xi_sens <- xi_sens[ok_sens]
  ni_sens <- ni_sens[ok_sens]
  lab_sens <- lab_sens[ok_sens]
  
  if (length(unique(tmp_sens$study_id[ok_sens])) >= 5L) {
    esc_sens <- escalc(measure = "PLO", xi = xi_sens, ni = ni_sens)
    esc_sens$study_id  <- factor(tmp_sens$study_id[ok_sens])
    esc_sens$effect_id <- seq_along(xi_sens)
    
    m3_sens <- rma.mv(yi, vi,
                      random = list(~1 | effect_id, ~1 | study_id),
                      data   = esc_sens, method = "REML")
    m2_sens <- rma.mv(yi, vi,
                      random = ~1 | study_id,
                      data   = esc_sens, method = "REML")
    
    prop_sens     <- inv_logit(m3_sens$b)
    propCI_sens   <- inv_logit(c(m3_sens$ci.lb, m3_sens$ci.ub))
    se_logit_sens <- sqrt(m3_sens$vb[1, 1])
    se_prop_sens  <- se_logit_sens * prop_sens * (1 - prop_sens)
    pred_sens     <- predict(m3_sens)
    pred_int_sens <- inv_logit(c(pred_sens$pi.lb, pred_sens$pi.ub))
    
    AIC_2_sens <- AIC(m2_sens); AIC_3_sens <- AIC(m3_sens)
    BIC_2_sens <- BIC(m2_sens); BIC_3_sens <- BIC(m3_sens)
    dAIC_sens  <- AIC_2_sens - AIC_3_sens
    dBIC_sens  <- BIC_2_sens - BIC_3_sens
    
    tau2_w_sens <- m3_sens$sigma2[1]; tau2_b_sens <- m3_sens$sigma2[2]
    n_eff_sens  <- nrow(esc_sens)
    inv_v_sens  <- 1/esc_sens$vi
    typ_sv_sens <- ((n_eff_sens - 1) * sum(inv_v_sens)) /
      (sum(inv_v_sens)^2 - sum(inv_v_sens^2))
    tot_v_sens  <- typ_sv_sens + tau2_w_sens + tau2_b_sens
    I2_w_sens   <- 100 * tau2_w_sens / tot_v_sens
    I2_b_sens   <- 100 * tau2_b_sens / tot_v_sens
    
    heterogeneity_log_sens <- rbind(heterogeneity_log_sens, data.frame(
      Threshold      = thr,
      Pooled         = prop_sens,
      SE_pooled      = se_prop_sens,
      CI_low         = propCI_sens[1],
      CI_high        = propCI_sens[2],
      PI_low         = pred_int_sens[1],
      PI_high        = pred_int_sens[2],
      Tau2_within    = tau2_w_sens,
      Tau2_between   = tau2_b_sens,
      I2_within      = I2_w_sens,
      I2_between     = I2_b_sens,
      AIC_2level     = AIC_2_sens,
      AIC_3level     = AIC_3_sens,
      dAIC           = dAIC_sens,
      BIC_2level     = BIC_2_sens,
      BIC_3level     = BIC_3_sens,
      dBIC           = dBIC_sens,
      stringsAsFactors = FALSE
    ))
    
    # Sensitivity forest plot PDF
    ticks_p <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.99)
    ticks_l <- log(ticks_p / (1 - ticks_p))
    pdf(paste0("Forest_", gsub("\\.", "_", thr), "_sens.pdf"), 9, 7)
    par(mar = c(5, 4, 5, 2))
    metafor::forest(m3_sens, slab = lab_sens,
                    atransf = inv_logit,
                    at      = ticks_l,
                    xlab    = "Prevalence of sub-optimal response (sensitivity)",
                    xlim    = c(min(ticks_l) - 1, max(ticks_l) + 1),
                    cex     = 0.8)
    axis(1, at = ticks_l,
         labels = formatC(ticks_p, format = "f", digits = 2), cex.axis = 0.8)
    title(main = paste("Sensitivity analysis – LDL goal", thr),
          cex.main = 1.1)
    legend("bottomright",
           legend = c(sprintf("I² within = %.1f%%",  I2_w_sens),
                      sprintf("I² between = %.1f%%", I2_b_sens),
                      sprintf("τ² total  = %.4f",    tau2_w_sens + tau2_b_sens),
                      sprintf("ΔAIC = %.2f",        dAIC_sens),
                      sprintf("ΔBIC = %.2f",        dBIC_sens),
                      sprintf("95%% PI %.3f–%.3f",  pred_int_sens[1], pred_int_sens[2])),
           bty = "n", cex = 0.8)
    dev.off()
    
    message("  ✓ Sensitivity analysis done (excluded zero numberofsors)")
  } else {
    message("  skipped sensitivity – fewer than 5 studies with numberofsors != 0")
  }
}

## After the loop finishes, save the sensitivity summary
write.csv(heterogeneity_log_sens, "heterogeneity_summary_sensitivity.csv", row.names = FALSE)
cat("\n✓ heterogeneity_summary_sensitivity.csv written\n")
print(heterogeneity_log_sens)

















###### sensitivity analysis
###############################################################################
# Sensitivity analysis – collapsed within-study effects (PLO method)
#   • Only thresholds defined by ≥5 studies
#   • Combine rows from same study within threshold
#   • Outputs: Pooled, SE, CI, PI, τ², I²
###############################################################################

library(metafor)

# --- 1. Load data ------------------------------------------------------------
df <- read.csv("C:\\Users\\ey270\\OneDrive - University of Cambridge\\Desktop\\Elias\\PHS\\Research Skills\\Dissertation Protocol\\data_3_7.csv",
               stringsAsFactors = FALSE, check.names = TRUE)
names(df) <- make.names(names(df))
df$Goal   <- as.character(df$Goal)

# --- 2. Extract LDL thresholds ------------------------------------------------
extract_threshold <- function(x) {
  m <- regmatches(x, regexpr("[0-9]+\\.?[0-9]*", x))
  if (length(m) == 0) return(NA_real_)
  as.numeric(m)
}
df$Goal_threshold <- vapply(df$Goal, extract_threshold, numeric(1))

# --- 3. Standardise study ID --------------------------------------------------
clean_study_id <- function(author) {
  a <- tolower(author)
  a <- gsub("[^a-z0-9]", "", a)
  a <- gsub("etal", "", a)
  a <- gsub("[0-9]{4}$", "", a)
  a
}
df$study_id <- vapply(df$author, clean_study_id, character(1))

# --- 4. Keep thresholds with ≥5 studies ---------------------------------------
df <- df[!is.na(df$Goal_threshold), ]
thr_counts <- table(unique(df[, c("Goal_threshold", "study_id")])$Goal_threshold)
eligible_thr <- as.numeric(names(thr_counts[thr_counts >= 5]))
if (!length(eligible_thr)) stop("No threshold has ≥5 studies")

cat("Thresholds retained:", paste(eligible_thr, collapse = ", "), "\n")

# --- 5. Helper functions ------------------------------------------------------
inv_logit <- function(x) exp(x) / (1 + exp(x))
fmt2 <- function(x) formatC(x, format = "f", digits = 2)
fmt3 <- function(x) formatC(x, format = "f", digits = 3)

summary_table <- data.frame()

# --- 6. Loop over eligible thresholds -----------------------------------------
for (thr in sort(eligible_thr)) {
  message("\n===== Analysing threshold ", thr, " mg/dL =====")
  tmp <- df[df$Goal_threshold == thr, ]
  
  # collapse within each study
  agg <- aggregate(cbind(numberofsors, samplesize) ~ study_id, data = tmp, sum, na.rm = TRUE)
  
  lab_lookup <- unique(tmp[, c("study_id", "author")])
  agg <- merge(agg, lab_lookup, by = "study_id", all.x = TRUE)
  agg$label <- agg$author
  
  xi <- as.numeric(agg$numberofsors)
  ni <- as.numeric(agg$samplesize)
  lab <- agg$label
  ok  <- !is.na(xi) & !is.na(ni) & xi >= 0 & ni > 0 & xi <= ni
  xi <- xi[ok]; ni <- ni[ok]; lab <- lab[ok]
  
  if (length(xi) < 5) {
    message("  skipped – fewer than 5 studies after collapsing")
    next
  }
  
  esc <- escalc(measure = "PLO", xi = xi, ni = ni)
  esc$study_id <- factor(agg$study_id[ok])
  
  m2 <- rma.mv(yi, vi, random = ~1 | study_id, data = esc, method = "REML")
  
  # pooled prevalence
  pooled     <- inv_logit(m2$b)
  pooled_CI  <- inv_logit(c(m2$ci.lb, m2$ci.ub))
  se_logit   <- sqrt(m2$vb[1, 1])
  se_pooled  <- se_logit * pooled * (1 - pooled)
  
  # prediction interval
  pred       <- predict(m2)
  pred_int   <- inv_logit(c(pred$pi.lb, pred$pi.ub))
  
  # heterogeneity stats
  tau2 <- m2$sigma2[1]
  inv_v  <- 1 / esc$vi
  typ_sv <- ((nrow(esc) - 1) * sum(inv_v)) / (sum(inv_v)^2 - sum(inv_v^2))
  I2    <- 100 * tau2 / (tau2 + typ_sv)
  
  cat(sprintf("  Pooled = %.3f (SE = %.3f) | 95%% CI = %.3f–%.3f | PI = %.3f–%.3f\n",
              pooled, se_pooled, pooled_CI[1], pooled_CI[2], pred_int[1], pred_int[2]))
  cat(sprintf("  τ² = %.4f | I² = %.1f%%\n", tau2, I2))
  
  # save row
  summary_table <- rbind(summary_table, data.frame(
    Threshold = thr,
    Pooled    = fmt3(pooled),
    SE        = fmt3(se_pooled),
    CI        = sprintf("(%s–%s)", fmt2(pooled_CI[1]), fmt2(pooled_CI[2])),
    PI        = sprintf("(%s–%s)", fmt2(pred_int[1]), fmt2(pred_int[2])),
    Tau2      = formatC(tau2, format = "f", digits = 4),
    I2        = fmt2(I2),
    stringsAsFactors = FALSE
  ))
  
  # forest plot
  ticks_p <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.99)
  ticks_l <- log(ticks_p / (1 - ticks_p))
  pdf(paste0("Forest_Sens_", thr, ".pdf"), 9, 7)
  par(mar = c(5,4,5,2))
  forest(m2, slab = lab,
         atransf = inv_logit,
         at = ticks_l,
         refline = NA,
         xlab = "Prevalence of sub-optimal response",
         xlim = c(min(ticks_l) - 1, max(ticks_l) + 1),
         cex = 0.8)
  axis(1, at = ticks_l,
       labels = formatC(ticks_p, format = "f", digits = 2), cex.axis = 0.8)
  title(main = paste("Sensitivity meta-analysis – LDL goal", thr), cex.main = 1.1)
  legend("bottomright",
         legend = c(
           sprintf("95%% PI = %s", sprintf("(%s–%s)", fmt2(pred_int[1]), fmt2(pred_int[2]))),
           sprintf("τ² = %.4f", tau2),
           sprintf("I² = %.1f%%", I2)),
         bty = "n", cex = 0.8)
  dev.off()
  cat(sprintf("  → Forest_Sens_%s.pdf saved\n", thr))
}

# --- 7. Save summary table ----------------------------------------------------
write.csv(summary_table, "sensitivity_summary.csv", row.names = FALSE)
cat("\n✓ sensitivity_summary.csv saved\n")
print(summary_table)







######################AGGREGATE SENSITIVITY ANALYSIS
###############################################################################
# Sensitivity analysis
#   • Combine all LDL‑C categories that belong to the same study
#   • One prevalence per study:  xi = Σ numberofsors,  ni = Σ samplesize
#   • Two‑level random‑effects (study as the only random effect)
#   • Logit‑transformed proportions (PLO)
###############################################################################
library(metafor)

## 1. Load data ---------------------------------------------------------------
df <- read.csv("C:\\Users\\ey270\\OneDrive - University of Cambridge\\Desktop\\Elias\\PHS\\Research Skills\\Dissertation Protocol\\data_3_7.csv",
               stringsAsFactors = FALSE, check.names = TRUE)

names(df) <- make.names(names(df))
df$Goal   <- as.character(df$Goal)

## 2. Build consistent study IDs (same helper you used) -----------------------
clean_study_id <- function(author) {
  a <- tolower(author)
  a <- gsub("[^a-z0-9]", "", a)
  a <- gsub("etal", "", a)
  a <- gsub("[0-9]{4}$", "", a)
  a
}
df$study_id <- vapply(df$author, clean_study_id, character(1))

## 3. Aggregate all LDL‑C categories inside each study ------------------------
agg <- aggregate(cbind(numberofsors, samplesize) ~ study_id,
                 data = df,
                 FUN = sum, na.rm = TRUE)

## 4. Add a label for the forest plot (author + year if available) ------------
# If your original data frame has "author" and "year", keep a lookup
authors <- unique(df[, c("study_id", "author")])
agg <- merge(agg, authors, by = "study_id", all.x = TRUE)

if ("year" %in% names(df)) {
  years <- unique(df[, c("study_id", "year")])
  agg <- merge(agg, years, by = "study_id", all.x = TRUE)
  agg$label <- with(agg, paste(author, year))
} else {
  agg$label <- agg$author
}

## 5. Validity filter ---------------------------------------------------------
xi <- as.numeric(agg$numberofsors)
ni <- as.numeric(agg$samplesize)
lab <- agg$label

ok <- !is.na(xi) & !is.na(ni) & xi >= 0 & ni > 0 & xi <= ni
xi  <- xi[ok];  ni <- ni[ok];  lab <- lab[ok]

if (length(unique(agg$study_id[ok])) < 5L) {
  stop("Fewer than 5 studies after aggregation – model not run")
}

## 6. Effect‑size calculation -------------------------------------------------
esc <- escalc(measure = "PLO", xi = xi, ni = ni)
esc$study_id <- factor(agg$study_id[ok])

## 7. Two‑level random‑effects model (only study level) -----------------------
m2 <- rma.mv(yi, vi, random = ~1 | study_id, data = esc, method = "REML")

## 8. Back‑transform pooled prevalence, CI, PI --------------------------------
inv_logit <- function(x) exp(x) / (1 + exp(x))

pooled   <- inv_logit(m2$b)
pooledCI <- inv_logit(c(m2$ci.lb, m2$ci.ub))
pred     <- predict(m2)
pred_int <- inv_logit(c(pred$pi.lb, pred$pi.ub))

cat(sprintf("\nPooled prevalence = %.3f  (95%% CI %.3f–%.3f)\n",
            pooled, pooledCI[1], pooledCI[2]))
cat(sprintf("95%% prediction interval = %.3f–%.3f\n\n",
            pred_int[1], pred_int[2]))

## 9. Forest plot -------------------------------------------------------------
ticks_p <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.99)
ticks_l <- log(ticks_p / (1 - ticks_p))      # logit ticks

pdf("Forest_Sensitivity_Aggregated.pdf", 9, 7)
par(mar = c(5,4,5,2))
forest(m2, slab = lab,
       atransf = inv_logit,
       at = ticks_l,
       refline = NA,
       xlab = "Prevalence of sub‑optimal response (aggregated)",
       xlim = c(min(ticks_l) - 1, max(ticks_l) + 1),
       cex = 0.8)
axis(1, at = ticks_l,
     labels = formatC(ticks_p, format = "f", digits = 2), cex.axis = 0.8)
title(main = "Sensitivity analysis – one prevalence per study", cex.main = 1.1)
legend("bottomright",
       legend = sprintf("95%% PI %.3f–%.3f", pred_int[1], pred_int[2]),
       bty = "n", cex = 0.8)
dev.off()
cat("→ Forest_Sensitivity_Aggregated.pdf saved\n")






















###############################################################################
# Sensitivity Analysis: Exclude studies by "Rafeeq" or "Kadikoylu"
###############################################################################

# Make a filtered copy of the data frame
df_sensitivity <- df[!grepl("Rafeeq|Kadikoylu", df$author, ignore.case = TRUE), ]

# Proceed with same analysis pipeline but use df_sensitivity
heterogeneity_log_sens <- data.frame()

for (thr in sort(unique(df_sensitivity$Goal_threshold))) {
  message("\n===== analysing threshold ", thr, " (sensitivity) =====")
  tmp <- df_sensitivity[df_sensitivity$Goal_threshold == thr, ]
  total_participants <- sum(as.numeric(tmp$samplesize), na.rm = TRUE)
  cat(sprintf("  total participants for %s mg/dL threshold: %d\n", thr, total_participants))
  
  xi <- suppressWarnings(as.numeric(tmp$numberofsors))
  ni <- suppressWarnings(as.numeric(tmp$samplesize))
  lab <- paste(tmp$author)
  lab2 <- paste(tmp$statins, tmp$Goal_details)
  
  ok <- !is.na(xi) & !is.na(ni) & xi >= 0 & ni > 0 & xi <= ni
  xi <- xi[ok]; ni <- ni[ok]; lab <- lab[ok]
  
  if (length(unique(tmp$study_id[ok])) < 5L) {
    message("  skipped – fewer than 5 studies with valid data")
    next
  }
  
  esc <- escalc(measure = "PLO", xi = xi, ni = ni)
  esc$study_id <- factor(tmp$study_id[ok])
  esc$effect_id <- seq_along(xi)
  
  m3 <- rma.mv(yi, vi,
               random = list(~1 | effect_id, ~1 | study_id),
               data = esc, method = "REML")
  m2 <- rma.mv(yi, vi,
               random = ~1 | study_id,
               data = esc, method = "REML")
  
  AIC_2 <- AIC(m2); AIC_3 <- AIC(m3)
  BIC_2 <- BIC(m2); BIC_3 <- BIC(m3)
  dAIC <- AIC_2 - AIC_3
  dBIC <- BIC_2 - BIC_3
  
  prop <- inv_logit(m3$b)
  propCI <- inv_logit(c(m3$ci.lb, m3$ci.ub))
  se_logit <- sqrt(m3$vb[1, 1])
  se_prop <- se_logit * prop * (1 - prop)
  pred <- predict(m3)
  pred_int <- inv_logit(c(pred$pi.lb, pred$pi.ub))
  
  cat(sprintf("  pooled prevalence %.3f (SE %.4f) 95%% CI %.3f–%.3f\n",
              prop, se_prop, propCI[1], propCI[2]))
  cat(sprintf("  prediction interval %.3f–%.3f\n", pred_int[1], pred_int[2]))
  cat(sprintf("  AIC 2‑level %.2f | 3‑level %.2f (Δ %.2f)  |  ",
              AIC_2, AIC_3, dAIC))
  cat(sprintf("BIC 2‑level %.2f | 3‑level %.2f (Δ %.2f)\n",
              BIC_2, BIC_3, dBIC))
  
  tau2_w <- m3$sigma2[1]; tau2_b <- m3$sigma2[2]
  n_eff <- nrow(esc)
  inv_v <- 1/esc$vi
  typ_sv <- ((n_eff - 1) * sum(inv_v)) /
    (sum(inv_v)^2 - sum(inv_v^2))
  tot_v <- typ_sv + tau2_w + tau2_b
  I2_w <- 100 * tau2_w / tot_v
  I2_b <- 100 * tau2_b / tot_v
  
  heterogeneity_log_sens <- rbind(heterogeneity_log_sens, data.frame(
    Threshold = thr,
    Pooled = prop,
    SE_pooled = se_prop,
    CI_low = propCI[1],
    CI_high = propCI[2],
    PI_low = pred_int[1],
    PI_high = pred_int[2],
    Tau2_within = tau2_w,
    Tau2_between = tau2_b,
    I2_within = I2_w,
    I2_between = I2_b,
    AIC_2level = AIC_2,
    AIC_3level = AIC_3,
    dAIC = dAIC,
    BIC_2level = BIC_2,
    BIC_3level = BIC_3,
    dBIC = dBIC,
    stringsAsFactors = FALSE
  ))
  
  pdf(paste0("Forest_Sensitivity_", gsub("\\.", "_", thr), ".pdf"), 9, 7)
  par(mar = c(5, 4, 5, 2))
  metafor::forest(m3, slab = lab, ilab=lab2,
                  atransf = inv_logit,
                  at = ticks_l,
                  xlab = "Prevalence of sub‑optimal response",
                  xlim = c(-13, 5),
                  cex = 0.8)
  axis(1, at = ticks_l,
       labels = formatC(ticks_p, format = "f", digits = 2), cex.axis = 0.8)
  title(main = paste("Sensitivity (excluded Rafeeq/Kadikoylu) – LDL goal", thr), cex.main = 1.1)
  
  dev.off()
  if (dev.cur() != 1) dev.off()
}

write.csv(heterogeneity_log_sens, "heterogeneity_summary_sensitivity.csv", row.names = FALSE)
cat("\n✓ heterogeneity_summary_sensitivity.csv written\n")
