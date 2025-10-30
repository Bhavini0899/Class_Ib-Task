gc()
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)
gse_data <- getGEO("GSE88837", GSEMatrix = TRUE)
expression_data <- exprs(gse_data$GSE88837_series_matrix.txt.gz)
feature_data <- fData(gse_data$GSE88837_series_matrix.txt.gz)
phenotype_data <- pData(gse_data$GSE88837_series_matrix.txt.gz)
sum(is.na(phenotype_data$source_name_ch1))
untar("Raw_Data/GSE88837_RAW.tar", exdir = "Raw_Data/CEL_Files")
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")
raw_data
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)
normalized_data <- rma(raw_data)
arrayQualityMetrics(expressionset = normalized_data[, c(1, 8, 9, 11, 20)],
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)
processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)
row_median <- rowMedians(as.matrix(processed_data))
row_median
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
threshold <- 1.8

abline(v = threshold, col = "blue", lwd = 1.5)
indx <- row_median > threshold
filtered_data <- processed_data[indx, ]
colnames(filtered_data) <- rownames(phenotype_data)
processed_data <- filtered_data
class(phenotype_data$source_name_ch1)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = "Visceral Adipose Tissue",
                 labels =  ("obessed"))
class(groups)
levels(groups)
