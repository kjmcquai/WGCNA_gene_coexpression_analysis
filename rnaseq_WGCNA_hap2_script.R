library(WGCNA)
library(tidyverse)

#------------------------------------------------------------------------------
# Bud Samples WGCNA
#------------------------------------------------------------------------------

# read in fpkm files for bud 
gene_fpkm_bud <- read.csv(file.choose(), row.names = 1)

# transpose table to match WGCNA format(samples as rows, genes as columns)
bud_t <- t(gene_fpkm_bud)
dim(bud_t)

# run goodSamples to make sure there are not too many missing genes
gsg <- goodSamplesGenes(bud_t, verbose = 3)
gsg$allOK

# run a sample tree to look for outliers
sampleTree <- hclust(dist(bud_t), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers (bud)", sub = "", xlab = "")

# determine which soft power to select (chart and plot)
powers <- c(1:20)
sft <- pickSoftThreshold(bud_t, powerVector = powers, verbose = 5)

plot(sft$fitIndices[,1], sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R^2)",
     type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.8, col = "red")

# run the program
bud_net <- blockwiseModules(bud_t,
                        power = 16,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = FALSE,
                        verbose = 3)

# view results
table(bud_net$colors)

pdf("bud_cluster_dendrogram_allblocks.pdf", width = 20, height = 8)
for (i in 1:length(bud_net$dendrograms)) {
  par(mar = c(2, 5, 3, 2))
  plotDendroAndColors(bud_net$dendrograms[[i]],
                      bud_net$colors[bud_net$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05,
                      main = paste("Block", i))
}
dev.off()


# add sex expression layer
#read in bud metadata
bud_meta <- read.csv(file.choose())

sex_bud <- as.data.frame(ifelse(bud_meta$sex == "F", 0, 1))
rownames(sex_bud) <- bud_meta$ids
colnames(sex_bud) <- "sex"
sex_bud

# fix name formatting to match metadata
rownames(bud_t) <- gsub("FPKM\\.", "", rownames(bud_t))
rownames(bud_t)

bud_MEs <- moduleEigengenes(bud_t, bud_net$colors)$eigengenes
bud_MEs <- orderMEs(bud_MEs)
bud_moduleTraitCor <- cor(bud_MEs, sex_bud, use = "p")
bud_moduleTraitPvalue <- corPvalueStudent(bud_moduleTraitCor, nrow(bud_t))


bud_textMatrix <- paste(signif(bud_moduleTraitCor, 2), "\n(",
                    signif(bud_moduleTraitPvalue, 1), ")", sep = "")
dim(bud_textMatrix) <- dim(bud_moduleTraitCor)

pdf("bud_module_trait_heatmap.pdf", width = 8, height = 14)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = bud_moduleTraitCor,
               xLabels = "Sex",
               yLabels = names(bud_MEs),
               ySymbols = names(bud_MEs),
               colorLabels = FALSE,
               colors = rev(blueWhiteRed(50)),
               textMatrix = bud_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-trait relationships\n(red = female-biased, blue = male-biased)")
dev.off()

# combine the heatmap and the dendrogram
pdf("bud_dendro_with_traits_allblocks.pdf", width = 20, height = 8)
for (i in 1:length(bud_net$dendrograms)) {
  geneModuleColors <- bud_net$colors[bud_net$blockGenes[[i]]]
  geneTraitCor <- bud_moduleTraitCor[paste0("ME", geneModuleColors), , drop = FALSE]
  geneTraitColors <- numbers2colors(geneTraitCor, signed = TRUE, colors = rev(blueWhiteRed(100)))
  plotColors <- cbind(geneModuleColors, geneTraitColors)
  colnames(plotColors) <- c("Module", "Sex correlation")
  
  par(mar = c(2, 5, 3, 2))
  plotDendroAndColors(bud_net$dendrograms[[i]],
                      plotColors,
                      groupLabels = colnames(plotColors),
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05,
                      main = paste("Block", i))
}
dev.off()

# determine significant modules 
sig_modules_bud <- data.frame(
  module = rownames(bud_moduleTraitCor),
  correlation = bud_moduleTraitCor[,1],
  pvalue = bud_moduleTraitPvalue[,1]
)
sig_modules_bud <- sig_modules_bud[order(sig_modules_bud$pvalue), ]
sig_modules_bud

write.csv(sig_modules_bud, "sig_modules_bud.csv")

# pull the gene IDs from the significant modules 
midnightblue_genes <- names(bud_net$colors)[bud_net$colors == "midnightblue"]
purple_genes <- names(bud_net$colors)[bud_net$colors == "purple"]
pink_genes <- names(bud_net$colors)[bud_net$colors == "pink"]
green_genes <- names(bud_net$colors)[bud_net$colors == "green"]
brown_genes <- names(bud_net$colors)[bud_net$colors == "brown"]

write.csv(midnightblue_genes, "bud_midnightblue_genes.csv", row.names = FALSE)
write.csv(purple_genes, "bud_purple_genes.csv", row.names = FALSE)
write.csv(pink_genes, "bud_pink_genes.csv", row.names = FALSE)
write.csv(green_genes, "bud_green_genes.csv", row.names = FALSE)
write.csv(brown_genes, "bud_brown_genes.csv", row.names = FALSE)

# filter to find which genes are present in the DEG list 
results_genes_bud <- read.csv(file.choose())

sig_degs_bud <- subset(results_genes_bud, qval <= 0.05)

sig_degs_bud$module <- bud_net$colors[match(sig_degs_bud$id, names(bud_net$colors))]

table(sig_degs_bud$module)

# save csvs of DEG genes in each module 
brown_degs <- subset(sig_degs_bud, module == "brown")
green_degs <- subset(sig_degs_bud, module == "green")
midnightblue_degs <- subset(sig_degs_bud, module == "midnightblue")
purple_degs <- subset(sig_degs_bud, module == "purple")
pink_degs <- subset(sig_degs_bud, module == "pink")

write.csv(brown_degs, "bud_brown_DEGs.csv", row.names = FALSE)
write.csv(green_degs, "bud_green_DEGs.csv", row.names = FALSE)
write.csv(midnightblue_degs, "bud_midnightblue_DEGs.csv", row.names = FALSE)
write.csv(purple_degs, "bud_purple_DEGs.csv", row.names = FALSE)
write.csv(pink_degs, "bud_pink_DEGs.csv", row.names = FALSE)

# Save the data to make cytoscape export easier
save(bud_net, bud_t, bud_MEs, bud_moduleTraitCor, bud_moduleTraitPvalue,
     sig_modules_bud, sex_bud,
     file = "bud_removeAF_wgcna_data.RData")

#------------------------------------------------------------------------------
# Flw Samples WGCNA
#------------------------------------------------------------------------------

# read in fpkm files for flw
gene_fpkm_flw <- read.csv(file.choose(), row.names = 1)

# transpose table to match WGCNA format(samples as rows, genes as columns)
flw_t <- t(gene_fpkm_flw)
dim(flw_t)

# run goodSamples to make sure there are not too many missing genes
gsg_flw <- goodSamplesGenes(flw_t, verbose = 3)
gsg_flw$allOK

# run a sample tree to look for outliers
sampleTree_flw <- hclust(dist(flw_t), method = "average")
plot(sampleTree_flw, main = "Sample clustering to detect outliers (flw)", sub = "", xlab = "")

# determine which soft power to select (chart and plot)
flw_powers <- c(1:20)
flw_sft <- pickSoftThreshold(flw_t, powerVector = flw_powers, verbose = 5)

plot(flw_sft$fitIndices[,1], flw_sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R^2)",
     type = "n")
text(flw_sft$fitIndices[,1], flw_sft$fitIndices[,2],
     labels = flw_powers, col = "red")
abline(h = 0.8, col = "red")

# run the program
flw_net <- blockwiseModules(flw_t,
                            power = 9,
                            TOMType = "unsigned",
                            minModuleSize = 30,
                            reassignThreshold = 0,
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            verbose = 3)

# view results
table(flw_net$colors)

pdf("flw_cluster_dendrogram_allblocks.pdf", width = 20, height = 8)
for (i in 1:length(flw_net$dendrograms)) {
  par(mar = c(2, 5, 3, 2))
  plotDendroAndColors(flw_net$dendrograms[[i]],
                      flw_net$colors[flw_net$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05,
                      main = paste("Block", i))
}
dev.off()

# add sex expression layer
#read in bud metadata
flw_meta <- read.csv(file.choose())

sex_flw <- as.data.frame(ifelse(flw_meta$sex == "F", 0, 1))
rownames(sex_flw) <- flw_meta$ids
colnames(sex_flw) <- "sex"
sex_flw

# fix name formatting to match metadata
rownames(flw_t) <- gsub("FPKM\\.", "", rownames(flw_t))
rownames(flw_t)

flw_MEs <- moduleEigengenes(flw_t, flw_net$colors)$eigengenes
flw_MEs <- orderMEs(flw_MEs)
flw_moduleTraitCor <- cor(flw_MEs, sex_flw, use = "p")
flw_moduleTraitPvalue <- corPvalueStudent(flw_moduleTraitCor, nrow(flw_t))


flw_textMatrix <- paste(signif(flw_moduleTraitCor, 2), "\n(",
                        signif(flw_moduleTraitPvalue, 1), ")", sep = "")
dim(flw_textMatrix) <- dim(flw_moduleTraitCor)

pdf("flw_module_trait_heatmap.pdf", width = 8, height = 14)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = flw_moduleTraitCor,
               xLabels = "Sex",
               yLabels = names(flw_MEs),
               ySymbols = names(flw_MEs),
               colorLabels = FALSE,
               colors = rev(blueWhiteRed(50)),
               textMatrix = flw_textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-trait relationships\n(red = female-biased, blue = male-biased)")
dev.off()

# combine the heatmap and the dendrogram
pdf("flw_dendro_with_traits_allblocks.pdf", width = 20, height = 8)
for (i in 1:length(flw_net$dendrograms)) {
  geneModuleColors <- flw_net$colors[flw_net$blockGenes[[i]]]
  geneTraitCor <- flw_moduleTraitCor[paste0("ME", geneModuleColors), , drop = FALSE]
  geneTraitColors <- numbers2colors(geneTraitCor, signed = TRUE, colors = rev(blueWhiteRed(100)))
  plotColors <- cbind(geneModuleColors, geneTraitColors)
  colnames(plotColors) <- c("Module", "Sex correlation")
  
  par(mar = c(2, 5, 3, 2))
  plotDendroAndColors(flw_net$dendrograms[[i]],
                      plotColors,
                      groupLabels = colnames(plotColors),
                      dendroLabels = FALSE,
                      hang = 0.03,
                      addGuide = TRUE,
                      guideHang = 0.05,
                      main = paste("Block", i))
}
dev.off()

# determine significant modules 
sig_modules_flw <- data.frame(
  module = rownames(flw_moduleTraitCor),
  correlation = flw_moduleTraitCor[,1],
  pvalue = flw_moduleTraitPvalue[,1]
)
sig_modules_flw <- sig_modules_flw[order(sig_modules_flw$pvalue), ]
sig_modules_flw

write.csv(sig_modules_flw, "sig_modules_flw.csv")

save(flw_net, flw_t, flw_MEs, flw_moduleTraitCor, flw_moduleTraitPvalue,
     sig_modules_flw, sex_flw,
     file = "flw_all_samples_wgcna_data.RData")