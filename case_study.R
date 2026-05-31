source("C:/GitHub/DepMap/distributional_knowledge_graph/fitting_functions.R")

# Assumes i have target and predictor loaded into memory
xp_data <- readRDS('C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/xp_filtered.rds')
chronos <- readRDS('C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/chronos_filtered.rds')
model <- read.csv('C:/GitHub/DepMap/data/26Q1/Model.csv')

# Describe each vectors first
x <- xp_data$PPARG
y <- chronos$RXRA_RXRB

# Case study
x_profile = summarize_vector_geometry(x,name='PPARG')
y_profile = summarize_vector_geometry(y,name='RXRA_RXRB')
xy_join = summarize_joint_geometry(x,y)

assoc <- summarize_global_linear_association(
  xp2$PPARG,
  x2$RXRA_RXRB
)

assoc

# Phase 3 asks: how does the expected value of Y change as X changes?
# This phase bridges cheap global summaries and full distributional models
# It answers if linearity is adequate
# Difference from phase 2:
# Not restricted to being linear
# explicitly directional
# starts modeling E[Y∣X]
# This is where Pearson correlation can fail
# What we'd like to detect here
# | relationship type | example                    |
#   | ----------------- | -------------------------- |
#   | linear            | steady decline             |
#   | saturating        | effect plateaus            |
#   | threshold         | only high PPARG matters    |
#   | U-shaped          | both extremes matter       |
#   | piecewise         | different slopes by regime |
#   | flat              | no mean structure          |

cond_mean_shape = summarize_conditional_mean_shape(x,y)

# Interpretation
# ΔR²: mostly linear <- Small Large ->  nonlinear structure matters
# ΔAIC: >0 - spline fit improved model, >>0 - strong evidence nonlinear structure exists
# nonlinearity_p: does spline explain variance beyond linear?
# monotonicity_score: 
# | score | meaning                    |
#   | ----- | -------------------------- |
#   | ~1    | strongly monotone          |
#   | ~0    | oscillatory / non-monotone |

# Phase 4 asks: how does variance of Y change as X changes?
# I.e. does X modulate uncertainty / dispersion / heterogeneity in Y?

linear_base = summarize_conditional_variance_structure(x,y, target_name = 'RXRA_RXRB', predictor_name = 'PPARG',mean_model = 'linear')
spline_base = summarize_conditional_variance_structure(x,y, target_name = 'RXRA_RXRB', predictor_name = 'PPARG',mean_model = 'spline', spline_df = 3)

# Phase 5
tail_summary <- summarize_tail_behavior(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB",
  left_threshold = -0.5
)

tail_summary

## Phase 6 Skewness
skew_summary <- summarize_skewness_asymmetry_structure(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB"
)

skew_summary

# Phase 7 - regime changes
regime_summary <- summarize_regime_threshold_structure(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB",
  left_tail_threshold = -0.5
)

regime_summary

# Visualize regime change
plot_regime_threshold(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  threshold = regime_summary$threshold,
  left_tail_threshold = -0.5,
  predictor_name = "PPARG expression",
  target_name = "RXRA_RXRB Chronos"
)

# Phase 8
shift_summary <- summarize_distributional_shift(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB"
)

shift_summary

# phase 9 - predictive utility
predictive_summary <- summarize_predictive_utility(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB",
  left_tail_threshold = -0.5
)

predictive_summary

# phase 10 - stability
stability_summary <- summarize_relationship_stability(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB",
  n_boot = 10,
  sample_frac = 0.80,
  left_tail_threshold = -0.5,
  seed = 1
)

stability_summary


subset(
  stability_summary,
  metric %in% c(
    "linear_slope",
    "left_tail_risk_difference",
    "regime_threshold",
    "regime_delta_aic",
    "wasserstein_1",
    "mean_shift"
  )
)

#########

pct_event_summary <- summarize_percentile_event_prediction(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB",
  event_prob = 0.20,
  direction = "left",
  n_folds = 5,
  spline_df = 3,
  seed = 1
)

pct_event_summary

top_event_summary <- summarize_percentile_event_prediction(
  predictor = xp2$PPARG,
  target = x2$RXRA_RXRB,
  predictor_name = "PPARG",
  target_name = "RXRA_RXRB",
  event_prob = 0.20,
  direction = "right"
)
