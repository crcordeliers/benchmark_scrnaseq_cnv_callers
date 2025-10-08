# fastCNV Integration Setup Guide

## Installation Steps

### 1. Create/Activate Conda Environment

```bash
# If environment doesn't exist yet, create it
micromamba env create -f conda_env_R.yml

# Activate the environment
micromamba activate r_env
```

### 2. Install fastCNV and Dependencies from GitHub/CRAN

Since Seurat (recent version) and fastCNV are not available on conda, install them directly in R:

```bash
R
```

Then in R:

```r
library(devtools)

# Install fastCNV and its data package
install_github("must-bioinfo/fastCNV")

# Seurat is needed by fastCNV (install latest version from CRAN)
# The conda version is very old, so install from CRAN
install.packages("Seurat")

# Verify installations
library(fastCNV)
library(Seurat)
library(data.table)

# Test that functions load
print("All packages loaded successfully!")
```

### 3. Verify Installation

Test that the fastCNV script can be parsed:

```bash
R --vanilla -e "parse('snakemake_pipeline/scripts/fastcnv.R'); print('fastcnv.R syntax OK')"
```

### 4. Test with a Dataset (if available)

Run fastCNV on a test dataset:

```bash
# Example: Run on SNU601 dataset
snakemake -s snakemake_pipeline/workflow.sm \
  results/output_SNU601/fastcnv/fastcnv_results.txt \
  --cores 1
```

## Troubleshooting

### Issue: Seurat version conflicts

If you encounter issues with Seurat, you may need to update it:

```r
remove.packages("Seurat")
install.packages("Seurat")
```

### Issue: Missing dependencies for fastCNV

If fastCNV installation fails, check the GitHub repo for required dependencies:
https://github.com/must-bioinfo/fastCNV

### Issue: CNV assay not found in Seurat object

This error means fastCNV did not create a CNV assay. Check:
1. That fastCNV ran successfully (check logs)
2. That the Seurat object has cells (not empty)
3. That reference cells are properly specified

## Testing Without Real Data

To create a minimal test dataset:

```r
# Create minimal test data
set.seed(123)
n_genes <- 1000
n_cells <- 50

# Create count matrix
counts <- matrix(rpois(n_genes * n_cells, lambda = 10),
                 nrow = n_genes, ncol = n_cells)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Cell", 1:n_cells)

# Create annotations
annotations <- data.frame(
  cell = colnames(counts),
  celltype = rep(c("cancer", "normal"), c(40, 10))
)

# Create ref groups
ref_groups <- data.frame(ref_groups = "normal")

# Save files
dir.create("data/input_test", recursive = TRUE)
write.table(cbind(gene = rownames(counts), counts),
            "data/input_test/count_matrix.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(annotations,
            "data/input_test/sample_annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(ref_groups,
            "data/input_test/ref_groups.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

Then run:

```bash
snakemake -s snakemake_pipeline/workflow.sm \
  results/output_test/fastcnv/fastcnv_results.txt \
  --cores 1
```

## Next Steps

Once fastCNV is running successfully:

1. Test on a real dataset (e.g., SNU601)
2. Run the full evaluation pipeline
3. Check that fastCNV appears in the evaluation outputs
4. Compare performance metrics with other methods

## File Locations

- **fastCNV wrapper script**: `snakemake_pipeline/scripts/fastcnv.R`
- **Reading functions**: `snakemake_pipeline/scripts/lib.R` (lines 506-581)
- **Snakemake rule**: `snakemake_pipeline/workflow.sm` (lines 11-28)
- **Evaluation integration**: `snakemake_pipeline/scripts/evaluate_results.R` (lines 119-129)
