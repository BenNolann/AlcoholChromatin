library(tximport)
library(DESeq2)
library(tximeta)
library(SummarizedExperiment)
library(readr)
library(dplyr)
library(GenomicFeatures)
library(AnnotationDbi)

##########################
### TXIMPORT
##########################

# Stage files for tximport
names <- c("UT_C", "UT_D", "UT_E", "E6_C", "E6_D", "E6_E", "E6R6_C","E6R6_D", "E6R6_E")
dir <- "/Zulu/bnolan/Projects/Personal/Ethanol/GEO/RNA-seq"
rep <- c(1,2,3)
files <- c(
  file.path(dir, paste0("Control_rep", rep, "SE_unstranded_anno_rsem.genes.results", sep = "")),
  file.path(dir, paste0("Ethanol_rep", rep, "SE_unstranded_anno_rsem.genes.results", sep = "")),
  file.path(dir, paste0("Ethanolwithdrawal_rep", rep, "SE_unstranded_anno_rsem.genes.results", sep = "")))

# Create metadata file for DESeq2
condition <- c(rep('UT',3), rep('E6',3), rep('E6R6',3))
replicate <- c(rep(c(1,2,3),3))

metadata <- data.frame(names, condition, replicate)
rownames(metadata) <- metadata$names

names(files) <- metadata$names


#prepare transcript id to gene symbol table
gtf = "~/Annotations/human/human/gtf/gencode/gencode.v47.annotation.gtf"
txdb = GenomicFeatures::makeTxDbFromGFF(gtf)


k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene$GENEID <- gsub("(\\.).*", "", tx2gene$GENEID) 

gtfTable <- "~/Annotations/human/human/gtf/gencode/gencode.v47_gene_annotation_table.txt"
gtfTabledf <- read.table(gtfTable, header = TRUE)


# Tximport
#txi <- tximport(files, type = "rsem", tx2gene = tx2gene, ignoreTxVersion = TRUE, )
#txi


# Use the gene ENSG* from the RSEM files directly
txi <- tximport(files, type = "rsem")
# Remove gene version/transcript
rownames(txi$counts) <- gsub("\\.\\d+$", "", rownames(txi$counts))
rownames(txi$abundance) <- gsub("\\.\\d+$", "", rownames(txi$abundance))
rownames(txi$length) <- gsub("\\.\\d+$", "", rownames(txi$length))


# 1. Sum Counts and Abundance using their own respective rownames
# This prevents dimension mismatch errors during sequential processing
txi$abundance <- rowsum(txi$abundance, rownames(txi$abundance))
txi$counts <- rowsum(txi$counts, rownames(txi$counts))

# 2. Get Max for Length
# Use aggregate to group by the rownames and apply max
len_df <- aggregate(txi$length, by = list(rownames(txi$length)), FUN = max)
rownames(len_df) <- len_df$Group.1
len_df$Group.1 <- NULL
txi$length <- as.matrix(len_df)

# 3. Build ENSG -> Symbol mapping
# Ensure ensg2symbol is pre-filtered to unique IDs
ensg2symbol <- gtfTabledf %>% 
  dplyr::select(Geneid, GeneSymbol) %>% 
  distinct() %>%
  mutate(Geneid = gsub("\\.\\d+$", "", Geneid)) %>%
  filter(!duplicated(Geneid))

# Create mapping vector
gene_map <- setNames(ensg2symbol$GeneSymbol, ensg2symbol$Geneid)

# 4. Create final names: prefer Symbol, fall back to ENSG if missing
raw_names <- ifelse(!is.na(gene_map[rownames(txi$counts)]), 
                    gene_map[rownames(txi$counts)], 
                    rownames(txi$counts))

# Ensure complete uniqueness for DESeq2
final_names <- make.unique(raw_names)

# 5. Apply final names to all matrices
rownames(txi$counts) <- final_names
rownames(txi$abundance) <- final_names
rownames(txi$length) <- final_names


keep_genes <- apply(txi$length, 1, min) > 0

# Subset the matrices to keep only the valid genes
txi$counts <- txi$counts[keep_genes, ]
txi$abundance <- txi$abundance[keep_genes, ]
txi$length <- txi$length[keep_genes, ]



# Verification
all(colnames(txi$counts) == rownames(metadata))



# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, metadata, design = ~ condition )


# Run DESeq2
dds <- DESeq(dds)

# Save object
setwd("/Zulu/bnolan/Projects/Personal/Ethanol/Assays/RNA-seq/ENCODE_CDE_4-21/")

saveRDS(dds, file="deseq2_dds.RDS")
#dds = readRDS("deseq2_dds.RDS")

##################
####### PCA (All genes)
##################
library(ggplot2)
rld <- rlog(dds)

data <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE)


percentVar <- round(100 * attr(data, "percentVar"))

data$replicate <- as.character(data$replicate)

ggplot(data, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  ggtitle("Ethanol treated HCT116 cells RNA-seq")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()



##################
####### Differential results
##################
res1 <- results(dds, contrast=c("condition","E6","UT"))
res2 <- results(dds, contrast=c("condition","E6R6","UT"))

summary(res1)
summary(res2)

nrow(res1[which(res1$log2FoldChange > 1 & res1$padj < 0.05),])
nrow(res1[which(res1$log2FoldChange < -1 & res1$padj < 0.05),])

nrow(res2[which(res2$log2FoldChange > 1 & res2$padj < 0.05),])
nrow(res2[which(res2$log2FoldChange < -1 & res2$padj < 0.05),])


##################
####### Volcano - Ethanol versus Control
##################

res1df <- as.data.frame(res1)

res1df <- res1df %>%
  dplyr::mutate(status = dplyr::case_when(
    abs(log2FoldChange) > log2(1.5) & padj < 0.05 ~ ">1.5-fold",
    
    TRUE ~ "NS"
  ))

ggplot(res1df, aes(
  x = log2FoldChange,
  y = -log10(padj),
  color = status
)) +
  geom_point() +
  scale_color_manual(
    values = c(
      ">1.5-fold" = "#004D40",
      "NS"        = "grey"
    ),
    name = "Fold Change"
  ) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 13)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_vline(
    xintercept = c(
      log2(1.5), -log2(1.5)
    ),
    linetype = "dotted"
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 20))
ggsave("volcano_E6.pdf", width = 6.22, height = 7) #4.18 x 4.6 in image


##################
####### Volcano - Ethanol Withdrawal versus Control
##################


res2df <- as.data.frame(res2)

res2df <- res2df %>%
  dplyr::mutate(status = dplyr::case_when(
    abs(log2FoldChange) > log2(1.5) & padj < 0.05 ~ ">1.5-fold",

    TRUE ~ "NS"
  ))

ggplot(res2df, aes(
  x = log2FoldChange,
  y = -log10(padj),
  color = status
)) +
  geom_point() +
  scale_color_manual(
    values = c(
      ">1.5-fold" = "#004D40",
      "NS"        = "grey"
    ),
    name = "Fold Change"
  ) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(0, 13)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_vline(
    xintercept = c(
      log2(1.5), -log2(1.5)
    ),
    linetype = "dotted"
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 20))
ggsave("volcano_E6R6.pdf", width = 6.22, height = 7) #4.18 x 4.6 in image


##################
####### Add metadata to results dataframes and merge Ethanol and Withdrawal
##################

gtfTabledf$Geneid <- gsub("\\.\\d+$", "", gtfTabledf$Geneid)
### Ethanol
match_symbol <- match(rownames(res1df), gtfTabledf$GeneSymbol)
match_id <- match(rownames(res1df), gtfTabledf$id)
final_idx <- ifelse(!is.na(match_symbol), match_symbol, match_id)
annotation_columns <- c("GeneSymbol", "Chromosome", "Start", "End", "Class", "Strand", "Length")
annotations_matched <- gtfTabledf[final_idx, annotation_columns]
res1df <- cbind(res1df, annotations_matched)

### Withdrawal
match_symbol <- match(rownames(res2df), gtfTabledf$GeneSymbol)
match_id <- match(rownames(res2df), gtfTabledf$id)
final_idx <- ifelse(!is.na(match_symbol), match_symbol, match_id)
annotation_columns <- c("GeneSymbol", "Chromosome", "Start", "End", "Class", "Strand", "Length")
annotations_matched <- gtfTabledf[final_idx, annotation_columns]
res2df <- cbind(res2df, annotations_matched)

# Merge Ethanol and Withdrawal together
res1df$gene_rowname <- rownames(res1df)
res2df$gene_rowname <- rownames(res2df)
annotation_cols <- c("GeneSymbol", "Chromosome", "Start", "End", "Class", "Strand", "Length")
join_keys <- c("gene_rowname", annotation_cols)
resall <- full_join(res1df, res2df, by = join_keys, suffix = c("_ethanol", "_withdrawal"))
rownames(resall) <- resall$gene_rowname
resall$gene_rowname <- NULL

# Remove genes with missing Data
getwd()
resallout <- resall %>% 
  na.omit(resall)

resallout %>%
  #dplyr::select(Chromosome, Start, End, Class, GeneSymbol.x, Strand, id, log2FoldChange.x, log2FoldChange.y, padj.x, padj.y) %>%
  write.table(file = "allnoNA_genes_16k.tsv", quote = FALSE, sep = "\t", row.names = TRUE)


##################
####### Create table for differential gene expression
##################


# Create a table for genes differential in both conditions and overlap of the >1.5-fold in 'status_withdrawal' and 'status_ethanol'
# Define log2 threshold corresponding to a 1.5 fold-change
# Change to lfc_thresh <- 1.5 if your criteria requires |log2FC| > 1.5 instead
# Define log2 threshold corresponding to a 1.5 fold-change
lfc_thresh <- log2(1.5)
total_genes <- nrow(resallout)

# Ethanol logical vectors
eth_up   <- !is.na(resallout$padj_ethanol) & resallout$padj_ethanol < 0.05 & resallout$log2FoldChange_ethanol > lfc_thresh
eth_down <- !is.na(resallout$padj_ethanol) & resallout$padj_ethanol < 0.05 & resallout$log2FoldChange_ethanol < -lfc_thresh
eth_tot  <- eth_up | eth_down

# Withdrawal logical vectors
with_up   <- !is.na(resallout$padj_withdrawal) & resallout$padj_withdrawal < 0.05 & resallout$log2FoldChange_withdrawal > lfc_thresh
with_down <- !is.na(resallout$padj_withdrawal) & resallout$padj_withdrawal < 0.05 & resallout$log2FoldChange_withdrawal < -lfc_thresh
with_tot  <- with_up | with_down

# Intersection logical vectors (significant in both conditions)
int_up   <- eth_up & with_up
int_down <- eth_down & with_down
int_tot  <- eth_tot & with_tot

# Count vectors
eth_counts  <- c(sum(eth_up), sum(eth_down), sum(eth_tot))
with_counts <- c(sum(with_up), sum(with_down), sum(with_tot))
int_counts  <- c(sum(int_up), sum(int_down), sum(int_tot))

# Main Summary Table (Counts + % of all dataset genes)
deg_summary <- data.frame(
  Condition      = c("Ethanol", "Withdrawal", "Intersect"),
  Up_Regulated   = c(eth_counts[1], with_counts[1], int_counts[1]),
  Down_Regulated = c(eth_counts[2], with_counts[2], int_counts[2]),
  Total_DEGs     = c(eth_counts[3], with_counts[3], int_counts[3]),
  Pct_Total_Genes = round(c(eth_counts[3], with_counts[3], int_counts[3]) / total_genes * 100, 2)
)

# Detailed Intersect Percentage Breakdown Table
intersect_pct_summary <- data.frame(
  Category               = c("Up-Regulated", "Down-Regulated", "Total DEGs"),
  Intersect_Count        = int_counts,
  Pct_of_Ethanol_DEGs    = round((int_counts / eth_counts) * 100, 2),
  Pct_of_Withdrawal_DEGs = round((int_counts / with_counts) * 100, 2),
  Pct_of_Total_Genes     = round((int_counts / total_genes) * 100, 2)
)

print(deg_summary)
print(intersect_pct_summary)

#write.table(resallout, file = "/Zulu/bnolan/Projects/Personal/Ethanol/Assays/RNA-seq/ENCODE_UT_E6_E6R6_all5samples_3-14-24/UsefulExpressionFiles/RNAseq_E6_E6R6_deseq2joinedtableNONA.tsv", quote = FALSE, row.names = FALSE)


##################
####### Pearson correlation of all genes with padj < 0.05
##################

# Subset genes to only include those with padj < 0.05. This is to visualize differential genes but some may want all padj values.
resalls <- resallout %>% filter(padj_ethanol < 0.05 | padj_withdrawal < 0.05)

resallslab = resalls %>%
  mutate(manE6 = ifelse(abs(log2FoldChange_ethanol) > log2(1.5), "diff", "notdiff"),
         manE6R6 = ifelse(abs(log2FoldChange_withdrawal) > log2(1.5), "diff", "notdiff"),
  ) 

resallslab <- resallslab %>%
  mutate(status = case_when(
    (manE6 == "notdiff" & manE6R6 == "notdiff") ~ "Not Affected",
    (manE6 == "notdiff" & manE6R6 == "diff") ~ "Delayed",
    (manE6 == "diff" & manE6R6 == "diff") ~ "Lasting",
    (manE6 == "diff" & manE6R6 == "notdiff") ~ "Acute"
  )) 

# Remove all NA padj genes.
resallslab = resallslab[complete.cases(resallslab),]

cor(resallslab$log2FoldChange_ethanol, resallslab$log2FoldChange_withdrawal)

################
# PCA of all padj < 0.05 genes across replicates (Control, Ethanol, Withdrawal)
################

### normalize counts
normalized_counts <- counts(dds, normalized = TRUE)
rownames(normalized_counts) <- gsub('\\..+$', '', rownames(normalized_counts))


normalized_counts <- as.data.frame(normalized_counts)
sig_genes <- rownames(normalized_counts) %in% rownames(resallslab)
normalized_sig_genes <- normalized_counts[sig_genes, ]

# Transpose the dataframe so columns become rows
df_transposed <- t(normalized_sig_genes)

# Standardize the data (recommended for PCA)
df_scaled <- scale(df_transposed)

# Perform PCA on transposed data
pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)

# View the summary of PCA results
summary(pca_result)

# Access PCA components
pca_result$rotation  # Loadings (eigenvectors)
pca_result$x         # Principal components (scores)

# Convert PCA results to a data frame
pca_df <- as.data.frame(pca_result$x)
pca_df$Variable <- colnames(normalized_sig_genes)  # Add original column names as labels

split_vars <- strsplit(pca_df$Variable, "_")

pca_df$condition <- sapply(split_vars, `[`, 1)
pca_df$replicate <- sapply(split_vars, `[`, 2)

ggplot(pca_df, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  ggtitle("nonNA padjusted genes (4781)")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
ggsave("PCA_noNApvalue_genes.pdf")

################
# Barplot of Lasting, Acute and Delayed genes
################

# Not affected, recovered, not recovered, worsened
viridis_colors <- c(
  "Not Affected" = "#440154",
  "Acute" = "#31688e",
  "Delayed" = "#fde725",
  "Lasting" = "#35b779"
)

# Barplot
resallslab %>%
  filter(status!="Not Affected") %>% 
  group_by(status) %>%
  summarise(n=n()) %>%
  #summarise(perc=n/sum(n)) %>%
  ggplot(mapping = aes(x=reorder(status, -n), y=n, fill=status, label=n)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = viridis_colors, limits = names(viridis_colors)) +
  geom_text(vjust="top", size = 10, color = "black", ) +
  ylab("") +
  xlab("") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        text = element_text(size = 30))
ggsave("gene_barplot.pdf", width = 10, height = 2) #4.18 x 4.6 in image


# Scatterplot
ggplot(resallslab, aes(x=log2FoldChange_ethanol, y=log2FoldChange_withdrawal, fill = factor(status))) +
  geom_point(shape=21, alpha = 0.6, size=3) +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
  xlab("log2FoldChange (EtOH / Control)") +
  ylab("log2FoldChange (With. / Control)") +
  geom_hline(yintercept = -log2(1.5), linetype="solid")+
  geom_hline(yintercept = log2(1.5), linetype="solid")+
  geom_vline(xintercept = -log2(1.5), linetype="solid")+
  geom_vline(xintercept = log2(1.5), linetype="solid")+
  scale_fill_manual(values = viridis_colors, limits = names(viridis_colors)) +
  theme_bw() +
  theme(
    legend.position="none",
    text = element_text(size = 15),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )


ggsave("gene_scatterplot_etoh_with.pdf", width = 6, height = 6)

################
# Table of differential genes
################

status_summary_table <- resallslab %>%
  group_by(manE6, manE6R6, status) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(pct_of_filtered = round((count / sum(count)) * 100, 2)) %>%
  rename(
    Ethanol_FC = manE6,
    Withdrawal_FC = manE6R6,
    Status = status,
    Count = count,
    Percent = pct_of_filtered
  ) %>%
  arrange(match(Status, c("Acute", "Delayed", "Lasting", "Not Affected")))

print(status_summary_table)

resallslab <- resalls %>%
  mutate(
    ethanol_change = case_when(
      !is.na(log2FoldChange_ethanol) & log2FoldChange_ethanol > log2(1.5)  ~ "gained",
      !is.na(log2FoldChange_ethanol) & log2FoldChange_ethanol < -log2(1.5) ~ "lost",
      TRUE ~ "unchanged"
    ),
    withdrawal_change = case_when(
      !is.na(log2FoldChange_withdrawal) & log2FoldChange_withdrawal > log2(1.5)  ~ "gained",
      !is.na(log2FoldChange_withdrawal) & log2FoldChange_withdrawal < -log2(1.5) ~ "lost",
      TRUE ~ "unchanged"
    ),
    status = case_when(
      ethanol_change != "unchanged" & withdrawal_change == "unchanged" ~ "Acute",
      ethanol_change == "unchanged" & withdrawal_change != "unchanged" ~ "Delayed",
      ethanol_change != "unchanged" & withdrawal_change != "unchanged" ~ "Lasting",
      ethanol_change == "unchanged" & withdrawal_change == "unchanged" ~ "Not Affected"
    )
  )

# 3. Detailed calculation table (shows directional combinations that form each status)
status_detail_table <- resallslab %>%
  group_by(status, ethanol_change, withdrawal_change) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percent = round((Count / sum(Count)) * 100, 2)) %>%
  rename(
    Status = status,
    Ethanol_Change = ethanol_change,
    Withdrawal_Change = withdrawal_change
  ) %>%
  arrange(match(Status, c("Acute", "Delayed", "Lasting", "Not Affected")))

# 4. Rolled-up table matching your exact status summary counts
status_summary_table <- status_detail_table %>%
  group_by(Status) %>%
  summarise(
    Total_Count = sum(Count),
    Total_Percent = sum(Percent),
    .groups = "drop"
  )

print(status_detail_table)
print(status_summary_table)


################
# Scatterplot with linear regression line and black/white
################

# 1. Define color/alpha logic
plot_data <- resallslab %>%
  mutate(
    PointColor = if_else(
      log2FoldChange_ethanol >= -0.58 & log2FoldChange_ethanol <= 0.58 &
        log2FoldChange_withdrawal >= -0.58 & log2FoldChange_withdrawal <= 0.58,
      "Center", "Outer"
    )
  )


# 1. Calculate the linear regression for the R2 value
fit <- lm(log2FoldChange_withdrawal ~ log2FoldChange_ethanol, data = resallslab)
r2_val <- round(summary(fit)$r.squared, 3)
label_text <- paste("R^2 == ", r2_val)

# 2. Generate the plot
ggplot(resallslab, aes(x = log2FoldChange_ethanol, y = log2FoldChange_withdrawal)) +
  # Blurred points
  geom_point(color = "black", alpha = 0.5, size = 1) +
  # Linear regression line
  geom_smooth(method = "lm", color = "grey", se = TRUE, linewidth = 1) +
  # Add dotted lines at 0.58 and -0.58 in all directions
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dotted", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = c(-0.58, 0.58), linetype = "dotted", color = "black", linewidth = 0.7) +
  # Add R2 value to the plot
  annotate("text", x = -1.5, y = 2.5, label = label_text, parse = TRUE, size = 6) +
  # Axes limits
  scale_x_continuous(limits = c(-2, 3)) +
  scale_y_continuous(limits = c(-2, 3)) +
  xlab("log2FoldChange (EtOH / Control)") +
  ylab("log2FoldChange (With. / Control)") +
  # Styling
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 15),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
ggsave('rnacorrelation.pdf')


#######################
# Output genes sorted by group (Acute / Lasting / Delayed ) then log fold change

rall <- resallslab %>%
  dplyr::select(Chromosome, Start, End, Class, GeneSymbol, Strand, log2FoldChange_ethanol, log2FoldChange_withdrawal, padj_ethanol, padj_withdrawal, status) %>%
  dplyr::arrange(status, -log2FoldChange_ethanol)
  #View()
  #write.table(file = "allsiggenes_status.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


##########
# Heatmap between E6 and E6R6 genes
##########
# Differentially expressed genes (Split across acute, lasting and delayed.)
resallslabAllgenes <- resallslab %>%
  dplyr::filter(
    (log2FoldChange_ethanol > log2(1.5) & padj_ethanol < 0.05) |
      (log2FoldChange_ethanol < -log2(1.5) & padj_ethanol < 0.05) |
      (log2FoldChange_withdrawal > log2(1.5) & padj_withdrawal < 0.05) |
      (log2FoldChange_withdrawal < -log2(1.5) & padj_withdrawal < 0.05)
  ) 


resallslabAllgenes %>%
  group_by(status) %>%
  summarize(n=n())

#1 Acute     947
#2 Delayed   792
#3 Lasting  1614


# Split these genes into E6, E6R6, Up in both, down in both, Alternating up and down.
resallslabAllgenes <- resallslabAllgenes %>%
  mutate(status = case_when(
    ((log2FoldChange_ethanol > log2(1.5)) & (log2FoldChange_withdrawal > log2(1.5))) ~ "BOTH_UP",
    ((log2FoldChange_ethanol < -log2(1.5)) & (log2FoldChange_withdrawal < -log2(1.5))) ~ "BOTH_DOWN",
    ((log2FoldChange_ethanol > log2(1.5)) & (log2FoldChange_withdrawal < log2(1.5) & log2FoldChange_withdrawal > -log2(1.5))) ~ "E6_UP",
    ((log2FoldChange_ethanol < -log2(1.5)) & (log2FoldChange_withdrawal < log2(1.5) & log2FoldChange_withdrawal > -log2(1.5))) ~ "E6_DOWN",
    ((log2FoldChange_withdrawal > log2(1.5)) & (log2FoldChange_ethanol < log2(1.5) & log2FoldChange_ethanol > -log2(1.5))) ~ "E6R6_UP",
    ((log2FoldChange_withdrawal < -log2(1.5)) & (log2FoldChange_ethanol < log2(1.5) & log2FoldChange_ethanol > -log2(1.5))) ~ "E6R6_DOWN",
  ))


resallslabAllgenes %>%
  group_by(status) %>%
  summarize(n=n())



# Create heatmap for ethanol and recovery genes

#create heatmap
library(pheatmap)

colors <- colorRampPalette(c("blue", "white", "red"))(100)


# Define colors: blue for above 0, red for below 0
# breaks <- seq(
#   max(abs(resallslabAllgenes$log2FoldChange.x), abs(resallslabAllgenes$log2FoldChange.y))*-log2(1.5), 
#   max(abs(resallslabAllgenes$log2FoldChange.x), abs(resallslabAllgenes$log2FoldChange.y)), 
#   length.out = length(colors))

breaks <- seq(-3,3,length.out = length(colors))



# Sort genes by E6 or E6R6 
resallslabAllgenes_sorted <- resallslabAllgenes[order(resallslabAllgenes$log2FoldChange_ethanol, decreasing=TRUE),]

# Sort genes by group (status)
resallslabAllgenes_sorted <- resallslabAllgenes %>%
  arrange(status)

# Subset and write sorted all genes to file
resallslabAllgenes_sorted_subset <-resallslabAllgenes_sorted %>%
  dplyr::select(
    GeneSymbol,
    log2FoldChange_ethanol,
    log2FoldChange_withdrawal,
    status
  )

# Add normalized read count from each sample to the dataframe
### normalize counts
normalized_counts <- counts(dds, normalized = TRUE)
rownames(normalized_counts) <- gsub('\\..+$', '', rownames(normalized_counts))
sig_genes <- rownames(normalized_counts) %in% resallslabAllgenes_sorted_subset$GeneSymbol
normalized_sig_genes <- normalized_counts[sig_genes, ]
normalized_counts_df <- as.data.frame(normalized_sig_genes)
normalized_counts_df_final <- normalized_counts_df[rownames(normalized_counts_df) %in%  resallslabAllgenes_sorted_subset$GeneSymbol,]
normalized_counts_df_final$GeneSymbol <- rownames(normalized_counts_df_final)

### Merge normalized counts with dataframe
# IF YOU WANT HEATMAP OF ALL REPLICATES USE THIS DF "resallslabAllgenes_sorted_subset_allcounts" TO DO IT
resallslabAllgenes_sorted_subset_allcounts <- merge(resallslabAllgenes_sorted_subset, normalized_counts_df_final, by = "GeneSymbol")


# Add columns for gene location from gtf file
#gtfTabledf$id <- gsub('\\..+$', '', gtfTabledf$id)
#gtfTabledf$Geneid <- NULL

#
resallslabAllgenes_sorted_subset_allcounts_locus <- merge(resallslabAllgenes_sorted_subset_allcounts, gtfTabledf, by = "GeneSymbol")

resallslabAllgenes_sorted_subset_allcounts_locus <- resallslabAllgenes_sorted_subset_allcounts_locus %>%
  dplyr::select(-GeneSymbol)


write.table(x=resallslabAllgenes_sorted_subset_allcounts_locus, file = "Allgenes_E6_E6recovery_counts.tsv", quote = FALSE, sep = "\t", row.names = FALSE)



resallslabAllgenes_sorted %>%
  group_by(status) %>%
  summarize(n=n()) %>%
  na.omit()

resallslabAllgenes_sorted_selected <- resallslabAllgenes_sorted %>%
  na.omit() %>%
  arrange(status, -log2FoldChange_ethanol) %>%
  dplyr::select(log2FoldChange_ethanol, log2FoldChange_withdrawal)

library(pheatmap)
# Create the heatmap
heat_out <- pheatmap(t(resallslabAllgenes_sorted_selected), 
                     color = colors, 
                     breaks = breaks,
                     show_rownames = FALSE, 
                     show_colnames = FALSE, 
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     scale = "none", 
                     legend = TRUE,
                     cluster_rows=FALSE,
                     cluster_cols=FALSE,
)     

output_file_path <- "heatmap_allgroups.pdf" 
pdf(output_file_path, height = 6.2, width=3.2)

heat_out
dev.off()



###########################
# GSEA REACTOME
###########################

#### ReactomePA
library(ReactomePA)

gsea_alcohol_raw <- gsePathway(geneList = gene_ranks_e6,
                               organism = "human",
                               minGSSize = 10,
                               maxGSSize = 500,
                               pvalueCutoff = 1,
                               pAdjustMethod = "BH")

dotplot_alcohol <- dotplot(gsea_alcohol_raw, showCategory = 20) + 
  ggtitle("Reactome GSEA: Alcohol")
dotplot_alcohol

gsea_withdrawal_raw <- gsePathway(geneList = gene_ranks_e6r6,
                                  organism = "human",
                                  minGSSize = 10,
                                  maxGSSize = 500,
                                  pvalueCutoff = 1,
                                  pAdjustMethod = "BH")

dotplot_withdrawal <- dotplot(gsea_withdrawal_raw, showCategory = 20) + 
  ggtitle("Reactome GSEA: Alcohol")
dotplot_withdrawal

# 1. Prepare data and filter for significance
df_alc <- as.data.frame(gsea_alcohol_raw)
df_with <- as.data.frame(gsea_withdrawal_raw)

# Keep only pathways where BOTH are p.adjust < 0.05
common_sig_paths <- intersect(
  df_alc$Description[df_alc$p.adjust < 0.05], 
  df_with$Description[df_with$p.adjust < 0.05]
)

# Subset and calculate mean NES for ordering
df_alc_sub <- df_alc %>% filter(Description %in% common_sig_paths) %>% select(Description, NES, p.adjust)
df_with_sub <- df_with %>% filter(Description %in% common_sig_paths) %>% select(Description, NES, p.adjust)

# Identify top 5 positive and negative based on Alcohol NES
top_pos <- df_alc_sub %>% filter(NES > 0) %>% slice_max(NES, n = 5)
top_neg <- df_alc_sub %>% filter(NES < 0) %>% slice_min(NES, n = 5)
final_paths <- c(top_pos$Description, top_neg$Description)

# Filter all data to these 10 and order by NES (Alcohol)
ordered_paths <- df_alc_sub %>% 
  filter(Description %in% final_paths) %>% 
  arrange(NES) %>% 
  pull(Description)

df_alc_sub <- df_alc_sub %>% filter(Description %in% final_paths)
df_with_sub <- df_with_sub %>% filter(Description %in% final_paths)

# Combine for plotting
plot_data_long <- rbind(
  df_alc_sub %>% mutate(Condition = "Alcohol"),
  df_with_sub %>% mutate(Condition = "Withdrawal")
) %>%
  mutate(Description = factor(Description, levels = ordered_paths))

# Prepare segment data
lollipop_segments <- df_alc_sub %>%
  rename(Alcohol = NES) %>%
  left_join(df_with_sub %>% rename(Withdrawal = NES), by = "Description") %>%
  mutate(Description = factor(Description, levels = ordered_paths))

# 2. Render the Lollipop Plot
ggplot() +
  geom_segment(
    data = lollipop_segments,
    aes(x = Alcohol, xend = Withdrawal, y = Description, yend = Description),
    color = "#cccccc", lwd = 1.0
  ) +
  geom_point(
    data = plot_data_long,
    aes(x = NES, y = Description, color = Condition),
    size = 4
  ) +
  scale_color_manual(
    values = c("Alcohol" = "#e41a1c", "Withdrawal" = "#377eb8"),
    name = "Condition"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#969696", lwd = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, margin = margin(t = 10)),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#f0f0f0"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  ) +
  labs(
    x = "Normalized Enrichment Score",
    title = "Top 10 Significant Reactome Pathways (Ordered by Alcohol NES)"
  )
ggsave('GSEA_alc_vwith_top10.pdf')



