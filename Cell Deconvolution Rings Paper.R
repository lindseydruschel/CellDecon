# Load required libraries
library(readxl)
library(dplyr)
library(biomaRt)
library(xCell)

## Next step, use xcell2 to use the MouseRNAseqData.xCell2Ref instead of human

# Step 1: Load the dataset
file_path <- "C:/Users/druschel/Downloads/Initial Dataset Rings Paper Transcriptomics.xlsx"
sheet_name <- "BioProbeCountMatrix"
data <- read_excel(file_path, sheet = sheet_name)

# Step 2: Define housekeeping genes and check their presence in the TargetName column
housekeeping_genes <- c("Actb", "Gapdh", "B2m", "Hprt", "Pgk1", "Ppia", "Rpl13a", "Tbp", "Ubc", "Ywhaz")
present_hk_genes <- data %>% filter(TargetName %in% housekeeping_genes)

# Step 3: If no housekeeping genes are found, stop
if (nrow(present_hk_genes) == 0) {
  stop("No housekeeping genes were found in the dataset.")
}

# Step 4: Extract sample columns and housekeeping gene data
sample_columns <- grep("ring", colnames(data), value = TRUE)  # Identify sample columns
hk_data <- present_hk_genes %>% select(TargetName, all_of(sample_columns))

# Step 5: Calculate geometric means for each sample and normalize the expression matrix
geo_means <- apply(hk_data[, -1], 2, function(x) exp(mean(log(x[x > 0]))))  # Avoid zeros
normalized_data <- data %>%
  mutate(across(all_of(sample_columns), ~ . / geo_means[match(cur_column(), names(geo_means))]))

# Step 6: Convert Mouse Genes to Human Orthologs using biomaRt
mouse_mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version = 110)
human_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 110)

# Extract mouse gene names from the dataset
mouse_genes <- normalized_data$TargetName

# Get human orthologs using biomaRt
mouse_to_human <- getBM(
  attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
  filters = "external_gene_name",
  values = mouse_genes,
  mart = mouse_mart
)

# Rename columns for clarity
colnames(mouse_to_human) <- c("Mouse_Gene", "Human_Gene")

# Filter out rows where no human ortholog is available
mouse_to_human <- mouse_to_human %>% filter(Human_Gene != "")

# Merge the human orthologs with the normalized dataset
normalized_data_human <- normalized_data %>%
  inner_join(mouse_to_human, by = c("TargetName" = "Mouse_Gene"))

# Step 7: Prepare for xCell Analysis
# Remove duplicates and set row names to human orthologs
expression_matrix_human <- normalized_data_human %>%
  group_by(Human_Gene) %>%
  summarise(across(where(is.numeric), mean)) %>%
  as.data.frame()  # Ensure it's a data.frame before setting rownames

# Set row names to "Human_Gene" and remove the column
rownames(expression_matrix_human) <- expression_matrix_human$Human_Gene
expression_matrix_human <- expression_matrix_human[, -1]  # Drop the "Human_Gene" column

# Ensure the matrix is numeric for xCell analysis
expression_matrix_human <- as.matrix(expression_matrix_human)

# Step 8: Run xCellAnalysis with the human orthologs
xcell_results <- xCellAnalysis(expression_matrix_human)

# Step 9: Output Results
print(head(xcell_results))

# Step 9: Visualize results with a heatmap
pheatmap(
  xcell_results,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "xCell Cell Type Enrichment Heatmap (FPKM)"
)

# OPTIONAL: Save results
write.csv(xcell_results, "C:/Users/druschel/Downloads/xCell_Results_Mouse_to_Human_Normalized.csv")
