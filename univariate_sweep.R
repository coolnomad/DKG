# Scan each individual feature first and assess its suitability for further analysis. 
# Drop features with weak distributions

source("C:/GitHub/DepMap/distributional_knowledge_graph/fitting_functions.R")

# Assumes i have target and predictor loaded into memory
xp_data <- readRDS('C:/GitHub/DepMap/data/26Q1/xp_26Q1.rds')
x <- readRDS('C:/GitHub/DepMap/data/26Q1/paralog_26Q1.rds')
model <- read.csv('C:/GitHub/DepMap/data/26Q1/Model.csv')


xp2 = xp_data[which(xp_data$ModelID %in% x$ModelID),]
xp2 = xp2[which(xp2$IsDefaultEntryForModel == 'Yes'),]
xp2 = xp2[which(xp2$IsDefaultEntryForMC == 'Yes'),]
x2 = x[which(x$ModelID %in% xp2$ModelID),]
xp2 = xp2[order(xp2$ModelID),]
table(xp2$ModelID == x2$ModelID)

colnames(xp2) <- sub("\\.\\..*$", "", colnames(xp2))

# Geometry filter to gene expression
xp_swept = sweep_vector_geometry(X = xp2[,-c(1:6)])

saveRDS(xp_swept, file='C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/xp_summary_stats.rds')

X_geometry_filtered <- filter_features_by_geometry(
  xp_swept,
  min_frac_complete = 0.95,
  min_n_unique = 20,
  min_frac_unique = 0.05,
  min_sd = 0.05,
  max_zero_frac = 0.95,
  min_effective_support_fraction = 0.40,
  max_bin_imbalance = 5
)

keepers <- c("ModelID",X_geometry_filtered$name[X_geometry_filtered$geometry_keep])

X_filtered <- xp2[, keepers, drop = FALSE]

table(X_geometry_filtered$geometry_keep)

sort(table(X_geometry_filtered$geometry_filter_reason), decreasing = TRUE)

head(
  X_geometry_filtered[order(X_geometry_filtered$sd, decreasing = TRUE), ],
  20
)

saveRDS(X_filtered,file='C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/xp_filtered.rds')

# Geometry filter to chronos response
chronos_swept = sweep_vector_geometry(X = x2[,-1])
saveRDS(chronos_swept, file='C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/chronos_summary_stats.rds')

chronos_geometry_filtered <- filter_features_by_geometry(
  chronos_swept,
  min_frac_complete = 0.95,
  min_n_unique = 20,
  min_frac_unique = 0.05,
  min_sd = 0.05,
  max_zero_frac = 0.95,
  min_effective_support_fraction = 0.40,
  max_bin_imbalance = 5
)

chronos_keepers <- c("ModelID",chronos_geometry_filtered$name[chronos_geometry_filtered$geometry_keep])

chronos_filtered <- x2[, chronos_keepers, drop = FALSE]

table(chronos_geometry_filtered$geometry_keep)

sort(table(chronos_geometry_filtered$geometry_filter_reason), decreasing = TRUE)

head(
  chronos_geometry_filtered[order(chronos_geometry_filtered$sd, decreasing = TRUE), ],
  20
)

saveRDS(chronos_filtered,file='C:/GitHub/DepMap/distributional_knowledge_graph/data/processed/chronos_filtered.rds')
