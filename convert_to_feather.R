data_dir <- "C:/GitHub/DepMap/data/26Q1"
out_dir  <- "C:/GitHub/DepMap/data/26Q1"

xp     <- readRDS(file.path(data_dir, "XP_26Q1.rds"))
crispr <- readRDS(file.path(data_dir, "CRISPR_26Q1.rds"))

cat("XP:    ", class(xp),     dim(xp),     "\n")
cat("CRISPR:", class(crispr), dim(crispr), "\n")
cat("XP rownames sample:    ", head(rownames(xp), 2),     "\n")
cat("CRISPR rownames sample:", head(rownames(crispr), 2), "\n")
cat("XP colnames sample:    ", head(colnames(xp), 3),     "\n")
cat("CRISPR colnames sample:", head(colnames(crispr), 3), "\n")

library(arrow)

# Convert matrix to data.frame with row names as first column
mat_to_df <- function(m) {
  df <- as.data.frame(m)
  df <- cbind(row_id = rownames(m), df)
  rownames(df) <- NULL
  df
}

cat("Writing XP feather...\n")
write_feather(mat_to_df(xp),     file.path(out_dir, "XP_26Q1.feather"))
cat("Writing CRISPR feather...\n")
write_feather(mat_to_df(crispr), file.path(out_dir, "CRISPR_26Q1.feather"))
cat("Done.\n")
