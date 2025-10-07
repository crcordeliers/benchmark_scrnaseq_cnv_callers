# ------------------------------------------------------------------------------
# Run fastCNV
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(fastCNV)
library(Seurat)
library(data.table)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_matrix <- snakemake@input$matrix
input_annotations <- snakemake@input$annot
input_ref_groups <- snakemake@input$ref_groups
input_gene_annotation<-snakemake@input$gene_pos

#Set the threshold for filtering expression dependent on 10X (0.1) vs Smartseq2 (1)
#Default the 10X threshold
param_expr_filter<-snakemake@params$expr_filter
if(is.null(param_expr_filter)){
  param_expr_filter<-0.1
}
param_expr_filter<-as.numeric(param_expr_filter)
  
output_file<-snakemake@output$cnv_file

output_file <- snakemake@output$cnv_file

threads <- as.numeric(snakemake@threads)
print(paste("Running fastCNV with", threads, "threads."))

# ------------------------------------------------------------------------------
print("Execute fastCNV")
# ------------------------------------------------------------------------------

# Create output directory in case it doesn't exist yet
output_dir <- dirname(output_file)
dir.create(output_dir, recursive=TRUE) # gives a warning if the directory exists already

# Load count matrix
data_matrix <- fread(input_matrix)
gene_names <- data_matrix$V1
data_matrix$V1 <- NULL
data_matrix <- as.matrix(data_matrix)
rownames(data_matrix) <- gene_names

# Load annotations
annotation <- fread(input_annotations, header=FALSE)
colnames(annotation) <- c("cell", "celltype")

# Read reference groups (saved in one tsv file)
ref_groups <- read.table(input_ref_groups, header=TRUE)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = data_matrix)

# Add annotations to metadata
cell_annotations <- annotation$celltype
names(cell_annotations) <- annotation$cell
seurat_obj[["celltype"]] <- cell_annotations

# Extract dataset name
dataset_name <- gsub("output_", "", unlist(strsplit(output_dir, split="/"))[2])

# Run fastCNV
seurat_obj <- fastCNV(seurat_obj,
                      sampleName = dataset_name,
                      referenceVar = "celltype",
                      referenceLabel = ref_groups$ref_groups,
                      printPlot = FALSE,
                      savePath = output_dir)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
