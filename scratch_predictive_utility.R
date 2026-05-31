summarize_predictive_utility <- function(
    predictor,
    target,
    predictor_name = "x",
    target_name = "y",
    n_folds = 5,
    spline_df = 3,
    left_tail_threshold = NULL,
    seed = 1
) {
  ok <- complete.cases(predictor, target)

  predictor <- predictor[ok]
  target <- target[ok]

  n <- length(predictor)

  if (n < 30 || length(unique(predictor)) < 10 || stats::sd(target) == 0) {
    return(data.frame(
      predictor = predictor_name,
      target = target_name,
      n = n,
      predictive_status = "insufficient_data",
      stringsAsFactors = FALSE
    ))
  }

  set.seed(seed)

  if (is.null(left_tail_threshold)) {
    left_tail_threshold <- unname(quantile(target, 0.10, na.rm = TRUE))
  }

  # ------------------------------------------------------
  # Fold assignment
  # ------------------------------------------------------

  fold_id <- sample(rep(seq_len(n_folds), length.out = n))

  linear_pred <- rep(NA_real_, n)
  spline_pred <- rep(NA_real_, n)

  tail_prob_linear    <- rep(NA_real_, n)
  tail_prob_linear_q20 <- rep(NA_real_, n)

  tail_event    <- target <= left_tail_threshold
  tail_event_q20 <- target <= unname(quantile(target, 0.20, na.rm = TRUE))

  # ------------------------------------------------------
  # Cross-validation loop
  # ------------------------------------------------------

  for (f in seq_len(n_folds)) {

    train_idx <- fold_id != f
    test_idx  <- fold_id == f

    x_train <- predictor[train_idx]
    y_train <- target[train_idx]

    x_test <- predictor[test_idx]

    # -----------------------------
    # Linear regression
    # -----------------------------

    linear_fit <- lm(y_train ~ x_train)

    linear_pred[test_idx] <- predict(
      linear_fit,
      newdata = data.frame(x_train = x_test)
    )

    # -----------------------------
    # Spline regression
    # -----------------------------

    spline_fit <- lm(
      y_train ~ splines::ns(x_train, df = spline_df)
    )

    spline_pred[test_idx] <- predict(
      spline_fit,
      newdata = data.frame(x_train = x_test)
    )

    # -----------------------------
    # Tail-event logistic model
    # -----------------------------

    glm_fit <- glm(
      tail_event[train_idx] ~ x_train,
      family = binomial()
    )

    tail_prob_linear[test_idx] <- predict(
      glm_fit,
      newdata = data.frame(x_train = x_test),
      type = "response"
    )

    glm_fit_q20 <- glm(
      tail_event_q20[train_idx] ~ x_train,
      family = binomial()
    )

    tail_prob_linear_q20[test_idx] <- predict(
      glm_fit_q20,
      newdata = data.frame(x_train = x_test),
      type = "response"
    )
  }

  # ------------------------------------------------------
  # Continuous prediction metrics
  # ------------------------------------------------------

  rmse <- function(obs, pred) {
    sqrt(mean((obs - pred)^2))
  }

  mae <- function(obs, pred) {
    mean(abs(obs - pred))
  }

  # null_rmse: RMSE of a mean-only (intercept-only) model — the baseline.
  # Provides the scale needed to interpret linear_rmse and spline_rmse.
  # Without null_rmse, a raw RMSE value has no reference: 0.3 could be
  # excellent or useless depending on the spread of Y.
  null_rmse <- stats::sd(target)

  linear_rmse <- rmse(target, linear_pred)
  spline_rmse <- rmse(target, spline_pred)

  linear_mae <- mae(target, linear_pred)
  spline_mae <- mae(target, spline_pred)

  linear_cor <- suppressWarnings(cor(target, linear_pred))
  spline_cor <- suppressWarnings(cor(target, spline_pred))

  linear_r2_cv <- linear_cor^2
  spline_r2_cv <- spline_cor^2

  # linear_skill_score: proportional reduction in RMSE vs the null model.
  # (null_rmse - linear_rmse) / null_rmse
  #  0 = no better than predicting the mean
  #  1 = perfect (zero RMSE)
  # Negative = worse than the null (possible under CV if the model overfits)
  # More interpretable than CV R² here because CV R² is correlation-based
  # (cor²) and does not penalize prediction bias; skill_score uses raw RMSE.
  linear_skill_score <- (null_rmse - linear_rmse) / null_rmse

  # ------------------------------------------------------
  # Tail-event prediction metrics
  # ------------------------------------------------------

  # AUROC
  calc_auc <- function(labels, scores) {

    ord <- order(scores, decreasing = TRUE)

    labels <- labels[ord]

    pos <- sum(labels)
    neg <- sum(!labels)

    if (pos == 0 || neg == 0) return(NA_real_)

    tpr <- cumsum(labels) / pos
    fpr <- cumsum(!labels) / neg

    sum(diff(c(0, fpr)) * (head(c(0, tpr), -1) + tail(c(0, tpr), -1)) / 2)
  }

  tail_auc <- calc_auc(tail_event, tail_prob_linear)

  # PR-AUC approximation
  calc_pr_auc <- function(labels, scores) {

    ord <- order(scores, decreasing = TRUE)

    labels <- labels[ord]

    tp <- cumsum(labels)
    fp <- cumsum(!labels)

    precision <- tp / (tp + fp + 1e-8)
    recall    <- tp / sum(labels)

    sum(diff(c(0, recall)) * tail(c(1, precision), -1))
  }

  tail_pr_auc <- calc_pr_auc(tail_event, tail_prob_linear)

  # Baseline prevalence
  tail_prevalence     <- mean(tail_event)
  tail_prevalence_q20 <- mean(tail_event_q20)

  # Calibration
  tail_brier     <- mean((tail_prob_linear     - tail_event)^2)
  tail_brier_q20 <- mean((tail_prob_linear_q20 - tail_event_q20)^2)

  tail_auc_q20    <- calc_auc(tail_event_q20,    tail_prob_linear_q20)
  tail_pr_auc_q20 <- calc_pr_auc(tail_event_q20, tail_prob_linear_q20)

  # pr_auc_lift: left_tail_pr_auc / left_tail_prevalence.
  # Normalizes PR-AUC against its random-classifier baseline, which equals
  # the event prevalence (not 0.5 as with AUROC). A lift of 2 means
  # precision-recall performance is 2× better than a random ranker.
  # Essential for comparing PR-AUC across pairs with different tail prevalences:
  # a PR-AUC of 0.25 is excellent when prevalence = 0.10 (lift = 2.5) but
  # mediocre when prevalence = 0.20 (lift = 1.25).
  pr_auc_lift     <- tail_pr_auc     / (tail_prevalence     + 1e-8)
  pr_auc_lift_q20 <- tail_pr_auc_q20 / (tail_prevalence_q20 + 1e-8)

  # ------------------------------------------------------
  # Output
  # ------------------------------------------------------

  data.frame(
    predictor = predictor_name,
    target    = target_name,
    n         = n,
    n_folds   = n_folds,

    null_rmse = null_rmse,

    linear_rmse = linear_rmse,
    spline_rmse = spline_rmse,
    delta_rmse  = linear_rmse - spline_rmse,

    linear_mae = linear_mae,
    spline_mae = spline_mae,
    delta_mae  = linear_mae - spline_mae,

    linear_cv_cor = linear_cor,
    spline_cv_cor = spline_cor,

    linear_cv_r2  = linear_r2_cv,
    spline_cv_r2  = spline_r2_cv,
    delta_cv_r2   = spline_r2_cv - linear_r2_cv,

    linear_skill_score = linear_skill_score,

    left_tail_threshold  = left_tail_threshold,
    left_tail_prevalence = tail_prevalence,

    left_tail_auc    = tail_auc,
    left_tail_pr_auc = tail_pr_auc,
    left_tail_brier  = tail_brier,
    pr_auc_lift      = pr_auc_lift,

    left_tail_prevalence_q20 = tail_prevalence_q20,
    left_tail_auc_q20        = tail_auc_q20,
    left_tail_pr_auc_q20     = tail_pr_auc_q20,
    left_tail_brier_q20      = tail_brier_q20,
    pr_auc_lift_q20          = pr_auc_lift_q20,

    predictive_status = "ok",
    stringsAsFactors = FALSE
  )
}

# Test
pred_util <- summarize_predictive_utility(
  x, y,
  predictor_name = "PPARG",
  target_name    = "RXRA_RXRB"
)
pred_util
