#!/usr/bin/env Rscript

# Load required libraries
library(CompoCor)
library(ggplot2)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    cat("Usage: Rscript cross_cor_SparXCC.r <OTU_TABLE> <EXPRESSION_MATRIX>\n")
    cat("Example: Rscript cross_cor_SparXCC.r merged_otu_table.tsv expression_matrix.tsv\n")
    quit(status = 1)
}

otu_file <- args[1]
expression_file <- args[2]

# Check if input files exist
if (!file.exists(otu_file)) {
    cat("Error: OTU table file not found:", otu_file, "\n")
    quit(status = 1)
}

if (!file.exists(expression_file)) {
    cat("Error: Expression matrix file not found:", expression_file, "\n")
    quit(status = 1)
}

cat("Reading input files...\n")
cat("OTU table:", otu_file, "\n")
cat("Expression matrix:", expression_file, "\n")

# READ DATA
raw_otu <- read.table(otu_file, header = TRUE)
raw_counts <- read.table(expression_file, header = TRUE)

cat("OTU table dimensions:", dim(raw_otu), "\n")
cat("Expression matrix dimensions:", dim(raw_counts), "\n")

# Split into day and night data
cat("Splitting data into day and night samples...\n")

counts_day <- raw_counts[, grep("day", names(raw_counts), value = TRUE)]
counts_night <- raw_counts[, grep("night", names(raw_counts), value = TRUE)]

otu_day_tmp <- raw_otu[, grep("day", names(raw_otu), value = TRUE)]
otu_night_tmp <- raw_otu[, grep("night", names(raw_otu), value = TRUE)]

# Fix colnames
colnames(counts_day) <- gsub("exp_day_", "", colnames(counts_day))
colnames(counts_night) <- gsub("exp_night_", "", colnames(counts_night))

colnames(otu_day_tmp) <- gsub("exp_day_", "", colnames(otu_day_tmp))
colnames(otu_night_tmp) <- gsub("exp_night_", "", colnames(otu_night_tmp))

cat("Day samples - OTU:", ncol(otu_day_tmp), "Expression:", ncol(counts_day), "\n")
cat("Night samples - OTU:", ncol(otu_night_tmp), "Expression:", ncol(counts_night), "\n")

# Sort column IDs to match between OTU and expression matrices
otu_day <- otu_day_tmp[, colnames(counts_day)]
otu_night <- otu_night_tmp[, colnames(counts_night)]

# Transpose matrices
t_otu_day <- t(otu_day)
t_otu_night <- t(otu_night)
t_counts_day <- t(counts_day)
t_counts_night <- t(counts_night)

# Check if sample IDs are identical between OTU and expression data
day_match <- identical(row.names(t_otu_day), row.names(t_counts_day))
night_match <- identical(row.names(t_otu_night), row.names(t_counts_night))

cat("Day sample IDs match:", day_match, "\n")
cat("Night sample IDs match:", night_match, "\n")

if (!day_match || !night_match) {
    cat("Warning: Sample IDs don't match between OTU and expression data\n")
}

# Find common samples between day and night
common_samples <- intersect(colnames(otu_day), colnames(otu_night))
cat("Common samples between day and night:", length(common_samples), "\n")

if (length(common_samples) == 0) {
    cat("Error: No common samples found between day and night data\n")
    quit(status = 1)
}

# Filter data to use only common samples
otu_day_common <- otu_day[, common_samples]
otu_night_common <- otu_night[, common_samples]
counts_day_common <- counts_day[, common_samples]
counts_night_common <- counts_night[, common_samples]

# Transpose filtered data
t_otu_day_common <- t(otu_day_common)
t_otu_night_common <- t(otu_night_common)
t_counts_day_common <- t(counts_day_common)
t_counts_night_common <- t(counts_night_common)

# Final check of column matching
day_common_match <- identical(colnames(t_otu_day_common), colnames(t_counts_day_common))
night_common_match <- identical(colnames(t_otu_night_common), colnames(t_counts_night_common))

cat("Day common samples column match:", day_common_match, "\n")
cat("Night common samples column match:", night_common_match, "\n")

# Run SparXCC for day samples
cat("Running SparXCC for day samples...\n")
sparxcc_day <- SparXCC_base(
    t_otu_day_common,
    t_counts_day_common,
    pseudo_count = 1,
    var_min = 1e-05,
    Find_m = TRUE,
    B_m = 100,
    cores = 8
)

# Save day results
write.table(sparxcc_day, file = "SparXCC_output_day_common_samples.txt", sep = "\t", row.names = T)
cat("Day results saved to: SparXCC_output_day_common_samples.txt\n")

# Run SparXCC for night samples
cat("Running SparXCC for night samples...\n")
sparxcc_night <- SparXCC_base(
    t_otu_night_common,
    t_counts_night_common,
    pseudo_count = 1,
    var_min = 1e-05,
    Find_m = TRUE,
    B_m = 100,
    cores = 8
)

# Save night results
write.table(sparxcc_night, file = "SparXCC_output_night_common_samples.txt", sep = "\t", row.names = T)
cat("Night results saved to: SparXCC_output_night_common_samples.txt\n")

# Create visualization if edge list files are available
create_visualizations <- function() {
    tryCatch({
        # Check if edge list files exist (assuming they are created by post-processing)
        day_edgelist_file <- "SparXCC_output_day_common_samples_edgelist.tsv"
        night_edgelist_file <- "SparXCC_output_night_common_samples_edgelist.tsv"
        
        if (file.exists(day_edgelist_file) && file.exists(night_edgelist_file)) {
            cat("Creating correlation distribution plots...\n")
            
            # Read edge list files
            sparxcc_day_edges <- read.csv(day_edgelist_file, sep = "\t", header = TRUE)
            sparxcc_night_edges <- read.csv(night_edgelist_file, sep = "\t", header = TRUE)
            
            # Add correlation type column
            sparxcc_day_edges$correlation_type <- ifelse(sparxcc_day_edges$CorrelationCoefficient > 0, "Positive", "Negative")
            sparxcc_night_edges$correlation_type <- ifelse(sparxcc_night_edges$CorrelationCoefficient > 0, "Positive", "Negative")
            
            # Set theme for plots
            theme_update(
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 14),
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
            )
            
            # Create day correlation distribution plot
            p_day <- ggplot(sparxcc_day_edges, aes(x = CorrelationCoefficient, fill = correlation_type)) +
                geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7) +
                scale_fill_manual(values = c("Positive" = "blue", "Negative" = "red")) +
                labs(
                    title = "Day",
                    x = "Correlation (SparXCC)",
                    y = "Frequency",
                    fill = "Correlation Sign"
                )
            
            ggsave("day_sparxcc_dist.png", plot = p_day, width = 10, height = 6, dpi = 300)
            
            # Create night correlation distribution plot
            p_night <- ggplot(sparxcc_night_edges, aes(x = CorrelationCoefficient, fill = correlation_type)) +
                geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7) +
                scale_fill_manual(values = c("Positive" = "blue", "Negative" = "red")) +
                labs(
                    title = "Night",
                    x = "Correlation (SparXCC)",
                    y = "Frequency",
                    fill = "Correlation Sign"
                )
            
            ggsave("night_sparxcc_dist.png", plot = p_night, width = 10, height = 6, dpi = 300)
            
            cat("Plots saved: day_sparxcc_dist.png, night_sparxcc_dist.png\n")
            
            # Print summary statistics
            cat("\n=== SUMMARY STATISTICS ===\n")
            cat("Day correlations:\n")
            cat("  Total:", nrow(sparxcc_day_edges), "\n")
            cat("  Positive:", sum(sparxcc_day_edges$correlation_type == "Positive"), "\n")
            cat("  Negative:", sum(sparxcc_day_edges$correlation_type == "Negative"), "\n")
            cat("  Mean correlation:", round(mean(sparxcc_day_edges$CorrelationCoefficient), 4), "\n")
            
            cat("\nNight correlations:\n")
            cat("  Total:", nrow(sparxcc_night_edges), "\n")
            cat("  Positive:", sum(sparxcc_night_edges$correlation_type == "Positive"), "\n")
            cat("  Negative:", sum(sparxcc_night_edges$correlation_type == "Negative"), "\n")
            cat("  Mean correlation:", round(mean(sparxcc_night_edges$CorrelationCoefficient), 4), "\n")
        } else {
            cat("Edge list files not found. Skipping visualization.\n")
            cat("To create plots, ensure edge list files are available:\n")
            cat("  -", day_edgelist_file, "\n")
            cat("  -", night_edgelist_file, "\n")
        }
    }, error = function(e) {
        cat("Error creating visualizations:", e$message, "\n")
    })
}

# Try to create visualizations
create_visualizations()

cat("\nSparXCC analysis completed successfully!\n")
cat("Output files:\n")
cat("  - SparXCC_output_day_common_samples.txt\n")
cat("  - SparXCC_output_night_common_samples.txt\n")
if (file.exists("day_sparxcc_dist.png")) {
    cat("  - day_sparxcc_dist.png\n")
    cat("  - night_sparxcc_dist.png\n")
}
