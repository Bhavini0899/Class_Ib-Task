data1 <- read.csv(file.choose())
data2 <- read.csv(file.choose())
write.csv(data1, "Raw Data/DEGs_Data_1.csv")
write.csv(data2, "Raw Data/DEGs_Data_2.csv")

input_dir <- "Raw Data"
if(!dir.exists(input_dir)){
  dir.create(input_dir)
}
output_dir <- "Results"
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")
result_list <- list()

classify_gene <- function(logFC, padj){
  ifelse(padj < 0.05 & logFC < -1, "Downregulated",
         ifelse(padj < 0.05 & logFC > 1, "Upregulated",
                "Non significant"))
}


for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
  if("padj" %in% names(data)){
    missing_count <- sum(is.na(data$padj))
    
    cat("Missing values in 'padj':", missing_count, "\n")
    data$padj[is.na(data$padj)] <- mean(data$padj, na.rm = TRUE)
  }
  
  
  if("logFC" %in% names(data)){
    missing_count <- sum(is.na(data$logFC))
    
    cat("Missing values in 'logFC':", missing_count, "\n")
    data$logFC[is.na(data$logFC)] <- mean(data$logFC, na.rm = TRUE)
  }
  
  data$Gene_Id <- classify_gene(data$logFC, data$padj)
  cat("Gene has been classified sucessfully !\n")
  result_list[[files]] <- data
  
  output_file_path <- file.path(output_dir, paste0("Classification"))
  write.csv(data, output_file_path, row.names = FALSE)
  
  gene_count <- table(data$Gene_Id)
  cat("summary counts for", file_names, "\n")
  print(gene_count)
}

result_1 = result_list[[1]]
result_2 = result_list[[2]]

save.image(file = "Bhavini_2_Assignment.Rdata")
