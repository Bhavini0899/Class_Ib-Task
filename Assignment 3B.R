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
gse_data <- getGEO("GSE69223", GSEMatrix = TRUE)
expression_data <- exprs(gse_data$GSE69223_series_matrix.txt.gz)
feature_data <- fData(gse_data$GSE69223_series_matrix.txt.gz)
phenotype_data <- pData(gse_data$GSE69223_series_matrix.txt.gz)
sum(is.na(phenotype_data$source_name_ch1))
untar("Raw_Data/GSE69223_RAW.tar", exdir = "Raw_Data/CEL_Files")
raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Result/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)
normalized_data <- rma(raw_data)
arrayQualityMetrics(expressionset = normalized_data[, c(1, 8, 9, 11, 20)],
                    outdir = "Result/QC_Normalized_Data",
                    force = TRUE)
processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)
row_median <- rowMedians(as.matrix(processed_data))
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
threshold <- 2.8
abline(v = threshold, col = "red", lwd = 2)
indx <- row_median > threshold
filtered_data <- processed_data[indx, ]
colnames(filtered_data) <- rownames(phenotype_data)
processed_data <- filtered_data
class(phenotype_data$source_name_ch1)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("PCa", "adj_normal"),
                 labels =  c("cancer", "normal"))
class(groups)
levels(groups)
save(feature_data, phenotype_data, processed_data, raw_data, file = "GSE69223.RData")
