# Loads transcripts abundances from salmon quant and aggregates values by genes
# Measures normalised TPM for downstream differential analyses.
# 20190604, cmdoret

library(tximport)
library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(apeglm)
library(tools)
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(biomaRt)

### PARSE CL ARGS ###

args <- commandArgs(trailingOnly=T)
# sample description file
samples <- read_tsv(args[1], comment='#')
# Sequencing units description file
units <- read_tsv(args[2], comment='#')
# Salmon output dir
salmon_dir <- args[3]
# Output file name
outfile <- args[4]

# Samples with RNAseq library
rna_samples <- units[units$libtype=='rnaseq', c('sample', 'replicate', 'batch')]
samples <- merge(samples, rna_samples, by='sample')
### LOAD DATA ###

quant_files <- file.path(salmon_dir, samples$sample, samples$replicate, "quant.sf")
# Include replicate number in sample names to prevent duplicated entries
samples$sample <- paste(samples$sample, samples$replicate, sep="_")
names(quant_files) <- samples$sample

#Get mapping between transcripts and genes
txdb <- as.data.frame(transcripts(EnsDb.Hsapiens.v86))
tx2gene <- dplyr::select(txdb, "tx_id", "gene_id")
rownames(tx2gene) <- NULL
colnames(tx2gene) <- c("TXNAME", "GENEID")

# Import transcripts counts from salmon into DESeq compatible object
txi.salmon <- tximport(
    quant_files, 
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance='lengthScaledTPM',
    ignoreTxVersion=T, 
)
rownames(samples) <- colnames(txi.salmon$counts)

# Differential expression analysis
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~ batch + condition)
# Pre-filtering genes with extremely low counts (<10 reads in total)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds, minReplicatesForReplace=Inf)

res <- results(dds)
# Computing log fold change with shrinkage estimator from apeglm
resLFCshrink <- lfcShrink(dds,
                    coef="condition_integration_vs_control",
                    type="apeglm")


### SAVE OUTPUT ###

base_out <- file_path_sans_ext(outfile)
# Writing table of gene counts for each sample
gene_counts <- data.frame(counts(dds))

# Convert ENSEMBL gene IDs to actual gene symbols
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", ensemblRedirect = FALSE)
symbols <- getBM(attributes=c('hgnc_symbol', "ensembl_gene_id"),
      filters="ensembl_gene_id",
      values=rownames(gene_counts),
      mart=ensembl)

# Add gene symbols to dataframes
add_symbols <- function(df, symbols){
    df <- df %>%
    rownames_to_column('ens_id') %>%
    left_join(symbols, by=c("ens_id" = "ensembl_gene_id")) %>%
    column_to_rownames("ens_id") %>%
    dplyr::select(hgnc_symbol, everything())
    return(df)
}
symbols$hgnc_symbol[symbols$hgnc_symbol == ""] <- symbols$ensembl_gene_id[symbols$hgnc_symbol == ""]
symbols <- symbols[!duplicated(symbols$ensembl_gene_id),]
gene_counts <- add_symbols(gene_counts, symbols)
resLFCshrink_df <- add_symbols(as.data.frame(resLFCshrink), symbols)

write.table(gene_counts,
            file=paste0(base_out, "_counts.tsv"),
            sep='\t', quote=F, row.names=T)

# Writing tablel of LFC and pvalues for all genes
write.table(resLFCshrink_df, file=outfile, sep='\t', quote=F, row.names=T)

### EXPLORATORY VISUALISATIONS ###

# Computing sample PCA
#vsd <- vst(dds, blind=FALSE)
# Saving PC1-PC2 plot
vsd <- vst(dds)
pca <- prcomp(t(assay(vsd)), scale=F)
svg(paste0(base_out, '_PCA_batch_uncorrected'))
plotPCA(vsd, 'batch')

dev.off()
# Correct batch effect
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
pca <- prcomp(t(assay(vsd)), scale=F)
#fviz_eig(pca)
svg(paste0(base_out, "_PCA_ind.svg"))
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
dev.off()
svg(paste0(base_out, "_PCA_var.svg"))
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             geom='arrow',     # Avoid text overlapping
             select.var=list(contrib=50),
             title="Top 50 genes contributing to PCs"
             )
dev.off()
# Computing pairwise sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$cell_type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
# Saving heatmap of sample distances
svg(paste0(base_out, "_dist_mat.svg"))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors ,main='Pairwise sample distance matrix')
dev.off()

# Saving volcano plot (LFC vs pvalue)
svg(paste0(base_out, "_volcano.svg"))
plot(res$log2FoldChange,
     -log10(res$padj),
     ylim=c(0, 20),
     main="Differential expression between integrated and control cancer cell lines",
     xlab="Log2 Fold change Integration / Control",
     ylab="-log10 adjusted P-Value")

dev.off()

# same with shrinked lfc
svg(paste0(base_out, "_volcano_LFCshrinked.svg"))
plot(resLFCshrink$log2FoldChange,
     -log10(resLFCshrink$padj),
     ylim=c(0, 20),
     main="Differential expression between integrated and control cancer cell lines",
     xlab="Log2 Fold change Integration/ Control",
     ylab="-log10 adjusted P-Value")

dev.off()

# Plotting Fold change vs coverage with raw LFC...
svg(paste0(base_out, "_plotMA.svg"))
plotMA(res)
dev.off()
# ...and with shrinked LFC
svg(paste0(base_out, "_plotMA_LFCshrink.svg"))
plotMA(resLFCshrink)
dev.off()

# Heatmap and clustering of gene counts for genes with top p-values
ordered_genes <- rownames(resLFCshrink)[order(resLFCshrink$padj)]
svg(paste0(base_out, "_heatmap_top300.svg"))
heatmap(as.matrix(gene_counts[ordered_genes[1:300],-1]),
        scale='row',
        labRow=gene_counts[ordered_genes[1:300],"hgnc_symbol"],
        main='Counts for 300 genes with top p-values')
dev.off()
