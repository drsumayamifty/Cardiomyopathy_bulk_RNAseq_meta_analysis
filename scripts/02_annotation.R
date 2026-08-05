library(genekitr)
library(tidyverse)
library(rio)

dir.create(
  "results/tables/annotated",
  showWarnings = FALSE,
  recursive = TRUE
)

geo_ids <- c("GSE71613", "GSE55296")

for (geo_id in geo_ids) {
  
  message("Processing ", geo_id, " ...")
  
  # Load DESeq2 results
  input_file <- paste0(
    "results/tables/DESeq2/",
    geo_id,
    ".csv"
  )
  
  if (!file.exists(input_file)) {
    message("Skipping ", geo_id, ": file not found")
    next
  }
  
  deseq_results <- rio::import(input_file)
  
  # Preserve original ID and remove Ensembl version suffix
  deseq_results <- deseq_results |>
    mutate(
      Gene_ID = as.character(Gene_ID),
      ENSEMBL = sub("\\..*$", "", Gene_ID)
    )
  
  # Annotate cleaned Ensembl IDs
  gene_info <- tryCatch(
    genekitr::genInfo(
      id = unique(deseq_results$ENSEMBL),
      org = "hs",
      unique = TRUE,
      keepNA = TRUE
    ),
    error = function(e) {
      message(
        "genInfo failed for ",
        geo_id,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(gene_info)) {
    next
  }
  
  # Join annotation without deleting DESeq2 rows
  annotated <- deseq_results |>
    left_join(
      gene_info,
      by = c("ENSEMBL" = "input_id")
    )
  
  # Rename annotation columns
  colnames(annotated)[colnames(annotated) == "symbol"] <- "Gene_Symbol" 
  colnames(annotated)[colnames(annotated) == "gene_name"] <- "Gene_Description" 
  colnames(annotated)[colnames(annotated) == "gene_biotype"] <- "Gene_Biotype"
  
  # Step 6: Keep only relevant columns 
  annotated <- annotated |> 
  select(any_of(c("Gene_ID", 
                  "ENSEMBL",
                  "Gene_Symbol", 
                  "Gene_Description", 
                  "Gene_Biotype", 
                  "baseMean", 
                  "log2FoldChange", 
                  "lfcSE", 
                  "stat", 
                  "pvalue", 
                  "padj"))) 
  # Step 7: Filter for protein-coding genes only 
  if ("Gene_Biotype" %in% colnames(annotated)) 
  { 
    annotated <- annotated |> 
      filter(Gene_Biotype == "protein_coding") 
  } 
  # Step 8: Export annotated results 
  output_file <- paste0("results/tables/annotated/", geo_id, ".csv") 
  export(annotated, output_file) 
  }