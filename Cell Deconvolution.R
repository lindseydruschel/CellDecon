# Load required libraries
library(xCell) # xCell2 allows custom refs; xCell uses built-in references
library(biomaRt)
library(readxl)
library(pheatmap)

#### On normal data ##### 

# Step 1: Load data
file_path <- "C:/Users/druschel/Downloads/Yulia Important Only.xlsx"
data <- read_excel(file_path, sheet = "Important Only")

# Step 2: Preprocess data (Keep gene names and expression columns)
data_clean <- data[, c("gene_name", grep("D2_No|D4_No", colnames(data), value = TRUE))]

# NEW Step 3: Convert Rat genes to Human orthologs using rat mart
# Connect to rat mart
mart <- useEnsembl(biomart = "ensembl", version = 110, dataset = "rnorvegicus_gene_ensembl")

# Query rat genes for Ensembl IDs and human orthologs
rat_orthologs <- getBM(
  attributes = c(
    "ensembl_gene_id",                     # Rat Ensembl Gene ID
    "external_gene_name",                  # Rat Gene Name
    "hsapiens_homolog_ensembl_gene",       # Human Ensembl Gene ID
    "hsapiens_homolog_associated_gene_name", # Human Gene Name
    "hsapiens_homolog_orthology_confidence" # Orthology Confidence
  ),
  filters = "external_gene_name",
  values = data_clean$gene_name,
  mart = mart
)

# Filter for high-confidence matches
rat_orthologs <- rat_orthologs[rat_orthologs$hsapiens_homolog_orthology_confidence == 1, ]

# Step 4: Merge with expression data
colnames(rat_orthologs) <- c("Rat_Ensembl_ID", "Rat_Gene", "Human_Ensembl_ID", "Human_Gene", "Orthology_Confidence")
data_human <- merge(data_clean, rat_orthologs, by.x = "gene_name", by.y = "Rat_Gene", all.x = TRUE)

# Step 5: Remove duplicates and handle missing values
data_human <- data_human[!duplicated(data_human$Human_Gene), ]  # Remove duplicates
data_human <- data_human[!is.na(data_human$Human_Gene), ]       # Remove rows with missing Human_Gene

# Assign row names to human orthologs
rownames(data_human) <- data_human$Human_Gene

# Step 6: Normalize data (Log2-transform, add +1 to avoid log(0) errors)
expression_matrix <- as.matrix(data_human[, grep("D2_No|D4_No", colnames(data_human))])
expression_matrix <- log2(expression_matrix + 1)

# Step 7: Run xCell analysis
xcell_results <- xCellAnalysis(expression_matrix)

# Step 8: Save results to CSV
write.csv(xcell_results, "C:/Users/druschel/Downloads/xCell_Results_Rat.csv", row.names = TRUE)

# Step 9: Visualize results with a heatmap
pheatmap(
  xcell_results,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "xCell Cell Type Enrichment Heatmap"
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
