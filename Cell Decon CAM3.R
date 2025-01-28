# Had to install git library (type this into console devtools::install_github("ChiungTingWu/CAM3"))
# NEEDS TO BE USING R 4.3.3 WITH BIOCONDUCTOR 3.19 AND A 64 BIT JAVA

# Load required libraries
library(CAM3)    # CAM3 package
library(readxl)  # For reading Excel files
library(ggplot2) # For visualization
library(dplyr)   # For data manipulation
library(gridExtra)

# Load your data from an Excel file
file_path <- "C:/Users/druschel/Downloads/Yulia all 17k important FPKM with MOUSE HOMOLOGs.xlsx"
bulk_data <- read_excel(file_path)

# Ensure gene IDs are unique
bulk_data$gene_id <- make.unique(as.character(bulk_data$gene_id))

# Convert tibble to data frame
bulk_data <- as.data.frame(bulk_data)

# Set row names and remove the `gene_id` column
rownames(bulk_data) <- bulk_data$gene_id
bulk_data <- bulk_data[, -which(colnames(bulk_data) == "gene_id")]

# Convert to matrix
data <- as.matrix(bulk_data)

# Step 1: Preprocess the data
dim.rdc <- 10  # Reduced dimension
thres.low <- 0.30  # Low threshold for filtering
thres.high <- 0.95  # High threshold for filtering
radius.thres <- 0.95  # Cosine radius for clustering
sim.thres <- 0.95  # Similarity threshold for merging clusters
cluster.num <- 50  # Minimum number of clusters
MG.num.thres <- 20  # Minimum number of marker genes

PrepResult <- CAM3Prep(
  data, dim.rdc, thres.low, thres.high,
  cluster.method = c("Fixed-Radius", "K-Means"),
  radius.thres, sim.thres, cluster.num, MG.num.thres
)

# Step 2: Marker gene clustering
fast.mode <- TRUE
MGResult <- CAM3MGCluster(PrepResult, fast.mode)

# Step 3: Estimate A and S matrices
generalNMF <- FALSE
ASestResult <- lapply(MGResult, CAM3ASest, PrepResult, data, generalNMF = generalNMF)

# Step 4: Run CAM3 integrated function
rCAM3 <- CAM3Run(data, K = 2:5, dim.rdc = 10, thres.low = 0.30, thres.high = 0.95)

# Output
print(rCAM3)

# Extract sample names from the column names of your bulk data matrix
sample_names <- colnames(bulk_data)

# Save the bar plot as a JPEG file
jpeg("cell_type_proportions.jpg", width = 800, height = 600)

# Create the bar plot with smaller text and adjusted margins
par(mar = c(8, 4, 4, 2))  # Adjust margins to fit rotated labels
barplot(
  t(A), 
  beside = TRUE, 
  col = rainbow(ncol(A)),
  legend.text = paste("Cell Type", 1:ncol(A)),
  args.legend = list(x = "topright", bty = "n"),
  main = "Predicted Cell Type Proportions",
  xlab = "Samples", 
  ylab = "Proportion",
  names.arg = sample_names,  # Use sample names as labels
  las = 2,  # Rotate labels for better visibility
  cex.names = 0.7  # Reduce font size of the sample names
)

# Close the graphics device
dev.off()

# Print confirmation message
cat("Bar plot saved as 'cell_type_proportions.jpg'\n")




##### Inspecting Marker Genes

# Load required library
library(openxlsx)  # Ensure you have the openxlsx package installed

# Inspect the marker gene matrix S
S <- ASestResult[[1]]@Sest  # Replace `[[1]]` with the relevant solution if needed

# Set a threshold to define "informative" genes (e.g., top genes with highest weights)
threshold <- 0.001  # Adjust based on your data
informative_genes <- apply(S, 2, function(x) rownames(S)[x > threshold])

# Assign names to the informative genes for clarity
names(informative_genes) <- paste("Cell Type", 1:ncol(S))

# Convert the list into a data frame with ranks
convert_to_ranked_df <- function(genes_list, top_n = NULL) {
  if (!is.null(top_n)) {
    # Limit the genes to top N for each cell type
    genes_list <- lapply(genes_list, function(genes) head(genes, top_n))
  }
  max_genes <- max(sapply(genes_list, length))  # Get the maximum number of genes across all cell types
  ranked_genes <- data.frame(
    Rank = 1:max_genes,
    do.call(cbind, lapply(genes_list, function(genes) {
      c(genes, rep(NA, max_genes - length(genes)))  # Fill shorter columns with NA
    }))
  )
  colnames(ranked_genes)[-1] <- names(genes_list)
  return(ranked_genes)
}

# Generate data frames for all markers, top 50, and top 15
all_markers_df <- convert_to_ranked_df(informative_genes)
top_50_df <- convert_to_ranked_df(informative_genes, top_n = 50)
top_15_df <- convert_to_ranked_df(informative_genes, top_n = 15)

# Write to an Excel file with three sheets
write.xlsx(
  list(
    "All Markers" = all_markers_df,
    "Top 50" = top_50_df,
    "Top 15" = top_15_df
  ),
  file = "Ranked_Informative_Genes.xlsx",
  rowNames = FALSE
)

cat("File saved as 'Ranked_Informative_Genes.xlsx' with sheets for all markers, top 50, and top 15.")
