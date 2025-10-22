gc()  #To clear memory

# checking annotation slot of the dataset
annotation(raw_data)
raw_data

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma", "AnnotationDbi", "affy", "hugene10sttranscriptcluster.db"))
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))

library(affy)
library(AnnotationDbi)
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(hugene10sttranscriptcluster.db)

cel_files <- list.files(pattern = "\\.cel", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
raw_data <- ReadAffy(filenames = cel_files)

cat("Platform:", annotation(raw_data), "\n")
cat("Samples:", length(sampleNames(raw_data)), "\n")


# Normalize with RMA
normalized_data <- rma(raw_data)
processed_data <- as.data.frame(exprs(normalized_data))
cat("Probes after normalization:", nrow(processed_data), "\n")

# Filter low-intensity probes
row_median <- rowMedians(as.matrix(processed_data))
threshold <- 4.0
keep_probes <- row_median > threshold
filtered_data <- processed_data[keep_probes, ]

cat("Probes after filtering:", nrow(filtered_data), "\n")

# Create groups (47 Control, 32 Alzheimer's)
sample_names <- sampleNames(raw_data)
groups <- c(rep("Control", 47), rep("Alzheimers_Disease", 32))

phenotype_data <- data.frame(
  Sample = sample_names,
  Group = groups
)


####Probe-to-Gene Mapping####
# making sure it is installed 
BiocManager::install("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)


probe_ids <- rownames(filtered_data)

gene_symbols <- mapIds(
  hugene10sttranscriptcluster.db,  
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

cat("Successfully mapped", sum(!is.na(gene_symbols)), "out of", length(probe_ids), "probes\n")


# Create gene mapping dataframe
gene_map_df <- data.frame(
  PROBEID = probe_ids,
  SYMBOL = gene_symbols
) %>% dplyr::filter(!is.na(SYMBOL))


# Check for multiple probes per gene
duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

cat("Multiple probes per gene summary:\n")
print(table(duplicate_summary$probes_per_gene))

# Merge annotation with expression data
processed_data_df <- filtered_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::left_join(gene_map_df, by = "PROBEID") %>%
  dplyr::filter(!is.na(SYMBOL))

# Collapse multiple probes per gene using average expression
expr_only <- processed_data_df %>% dplyr::select(-PROBEID, -SYMBOL)
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

cat("Final gene count:", nrow(averaged_data), "\n")

# Convert to matrix for limma
data <- as.matrix(averaged_data)


####Differential Expression Analysis####

groups <- factor(phenotype_data$Group,
                 levels = c("Control", "Alzheimers_Disease"),
                 labels = c("normal", "alzheimer"))

# Create design matrix
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# Fit linear model
fit_1 <- lmFit(data, design)

# Define contrast: Alzheimer's vs Control
contrast_matrix <- makeContrasts(alzheimer_vs_normal = alzheimer - normal,
                                 levels = design)

# Apply contrasts and compute statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)
fit_2 <- eBayes(fit_contrast)

# Extract differentially expressed genes
deg_results <- topTable(fit_2,
                        coef = "alzheimer_vs_normal",
                        number = Inf,
                        adjust.method = "BH")

# Classify DEGs (adj.P.Val < 0.05 and |logFC| > 1)
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))


# Subset DEGs
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")
deg_updown <- rbind(upregulated, downregulated)

cat("DEG SUMMARY:\n")
cat("- Upregulated genes:", nrow(upregulated), "\n")
cat("- Downregulated genes:", nrow(downregulated), "\n")
cat("- Total significant DEGs:", nrow(deg_updown), "\n")

if(!dir.exists("Results")) dir.create("Results")
if(!dir.exists("Result_Plots")) dir.create("Result_Plots")

# Save results
write.csv(deg_results, "Results/DEGs_Results_Alzheimer.csv")
write.csv(upregulated, "Results/Upregulated_DEGs_Alzheimer.csv")
write.csv(downregulated, "Results/Downregulated_DEGs_Alzheimer.csv")

# Volcano Plot
png("Result_Plots/volcano_plot_alzheimer.png", width = 2000, height = 1500, res = 300)
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No" = "grey")) +
  theme_minimal() +
  labs(title = "Alzheimer's Disease vs Control", x = "log2 Fold Change", y = "-log10(Adjusted P-value)")
dev.off()

# Heatmap of Top 25 DEGs
deg_updown_ordered <- deg_updown[order(deg_updown$adj.P.Val), ]
top_25_genes <- rownames(deg_updown_ordered)[1:min(25, nrow(deg_updown_ordered))]

genes_in_both <- intersect(top_25_genes, rownames(data))

if(length(genes_in_both) > 0) {
  heatmap_data <- data[genes_in_both, ] }
  
  png("Result_Plots/heatmap_top25_DEGs_alzheimer.png", width = 2000, height = 1500, res = 300)
  pheatmap(
    heatmap_data,
    scale = "row",
    main = paste("Top", length(genes_in_both), "Differentially Expressed Genes"),
    color = colorRampPalette(c("blue", "white", "red"))(100)
  )
  dev.off()

  
  
  cat("EXACT NUMBERS:\n")
  cat("Probes before collapsing:", nrow(filtered_data), "\n")
  cat("Genes after collapsing:", nrow(averaged_data), "\n")
  cat("Upregulated genes:", nrow(upregulated), "\n")
  cat("Downregulated genes:", nrow(downregulated), "\n")
  cat("Total DEGs:", nrow(deg_updown), "\n")
  
  
  