# Load required libraries
library(xCell2)   # xCell2 for custom references and built-in mouse references
library(biomaRt)  # biomaRt for querying orthologs
library(readxl)   # readxl for loading data
library(pheatmap) # pheatmap for heatmap visualization

#### On normal data ##### 

# Step 1: Load data
file_path <- "C:/Users/druschel/Downloads/Yulia Important Only.xlsx"
data <- read_excel(file_path, sheet = "Important Only")

# Step 2: Preprocess data (Keep gene names and expression columns)
data_clean <- data[, c("gene_name", grep("D2_No|D4_No", colnames(data), value = TRUE))]

# Step 3: Convert Rat genes to Mouse orthologs using rat mart
# Connect to rat mart
mart_rat <- useEnsembl(biomart = "ensembl", version = 110, dataset = "rnorvegicus_gene_ensembl")

# Query rat genes for Ensembl IDs and mouse orthologs
rat_orthologs <- getBM(
  attributes = c(
    "ensembl_gene_id",                     # Rat Ensembl Gene ID
    "external_gene_name",                  # Rat Gene Name
    "mmusculus_homolog_ensembl_gene",     # Mouse Ensembl Gene ID
    "mmusculus_homolog_associated_gene_name", # Mouse Gene Name
    "mmusculus_homolog_orthology_confidence" # Orthology Confidence
  ),
  filters = "external_gene_name",
  values = data_clean$gene_name,
  mart = mart_rat
)

# Filter for high-confidence matches
rat_orthologs <- rat_orthologs[rat_orthologs$mmusculus_homolog_orthology_confidence == 1, ]

# Step 4: Merge with expression data
colnames(rat_orthologs) <- c("Rat_Ensembl_ID", "Rat_Gene", "Mouse_Ensembl_ID", "Mouse_Gene", "Orthology_Confidence")
data_mouse <- merge(data_clean, rat_orthologs, by.x = "gene_name", by.y = "Rat_Gene", all.x = TRUE)

# Step 5: Remove duplicates and handle missing values
data_mouse <- data_mouse[!duplicated(data_mouse$Mouse_Gene), ]  # Remove duplicates
data_mouse <- data_mouse[!is.na(data_mouse$Mouse_Gene), ]       # Remove rows with missing Mouse_Gene

# Assign row names to mouse orthologs
rownames(data_mouse) <- data_mouse$Mouse_Gene

# Step 6: Normalize data (Log2-transform, add +1 to avoid log(0) errors)
expression_matrix <- as.matrix(data_mouse[, grep("D2_No|D4_No", colnames(data_mouse))])
expression_matrix <- log2(expression_matrix + 1)

# Step 7: Run xCell2 analysis using the built-in mouse references
xcell_results_fpkm_mouse <- xCell(expression_matrix)

# Step 8: Save results to CSV
write.csv(xcell_results_fpkm_mouse, "C:/Users/druschel/Downloads/xCell2_Results_FPKM_Mouse.csv", row.names = TRUE)

# Step 9: Visualize results with a heatmap
pheatmap(
  xcell_results_fpkm_mouse,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "xCell2 Mouse Cell Type Enrichment Heatmap (FPKM)"
)
#### Same thing but for FPKM Vals #####
# Load required libraries
library(xCell)
library(biomaRt)
library(readxl)
library(pheatmap)

# Step 1: Load data
file_path <- "C:/Users/druschel/Downloads/D2vsD4_deg w marker check (4).xlsx"
data <- read_excel(file_path, sheet = "D2vsD4_deg.xls")

# Step 2: Preprocess data (Keep gene names and FPKM columns)
data_fpkm <- data[, c("gene_name", grep("_fpkm$", colnames(data), value = TRUE))]

# Step 3: Convert Rat genes to Human orthologs using rat mart
mart <- useEnsembl(biomart = "ensembl", version = 110, dataset = "rnorvegicus_gene_ensembl")
rat_orthologs <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_confidence"
  ),
  filters = "external_gene_name",
  values = data_fpkm$gene_name,
  mart = mart
)

# Filter for high-confidence matches
rat_orthologs <- rat_orthologs[rat_orthologs$hsapiens_homolog_orthology_confidence == 1, ]

# Step 4: Merge with FPKM data
colnames(rat_orthologs) <- c("Rat_Ensembl_ID", "Rat_Gene", "Human_Ensembl_ID", "Human_Gene", "Orthology_Confidence")
data_human_fpkm <- merge(data_fpkm, rat_orthologs, by.x = "gene_name", by.y = "Rat_Gene", all.x = TRUE)

# Step 5: Remove duplicates and handle missing values
data_human_fpkm <- data_human_fpkm[!duplicated(data_human_fpkm$Human_Gene), ]
data_human_fpkm <- data_human_fpkm[!is.na(data_human_fpkm$Human_Gene), ]
rownames(data_human_fpkm) <- data_human_fpkm$Human_Gene

# Step 6: Use FPKM values directly
fpkm_columns <- grep("_fpkm$", colnames(data_human_fpkm))  # Identify FPKM columns
fpkm_matrix <- as.matrix(data_human_fpkm[, fpkm_columns])  # Subset FPKM values

# Step 7: Run xCell analysis
xcell_results_fpkm <- xCellAnalysis(fpkm_matrix)

# Step 8: Save results to CSV
write.csv(xcell_results_fpkm, "C:/Users/druschel/Downloads/xCell_Results_FPKM_D2vsD4.csv", row.names = TRUE)

# Step 9: Visualize results with a heatmap
pheatmap(
  xcell_results_fpkm,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "xCell Cell Type Enrichment Heatmap (FPKM)"
)
