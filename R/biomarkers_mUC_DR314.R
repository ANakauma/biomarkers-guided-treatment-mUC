# Author:   J. Alberto Nakauma-Gonzalez (j.nakaumagonzalez@erasmusmc.nl)
# Date:     11-02-2025
# Function: Script for the paper: Novel molecular biomarkers to guide treatment-decision making
# in metastatic urothelial cancer - A patient cohort analysis


# Load libraries ----------------------------------------------------------
pacman::p_load('plyr', 'dplyr', 'tidyr', 'ggplot2', 'pbapply', 'GSVA', 'reshape2')
library(consensusMIBC)

# path to data and output dir
path.hmf <- "/studycohorts/DR314/"
path.RData <- "/studycohorts/DR314/postHMF/RData/"
odir <- "/Projects/DR314_HMF/output/"
path.rna <- "/studycohorts/DR314/postHMF/RNAseq/counts_kallisto/"

# load pre-processed data (mutations and patient metadata)
load(paste0(path.RData,"DR314.MetaData.RData"))
load(paste0(path.RData,"results.Cohort.RData"))


# define colors for variables
color_TSEscore <-  c('TSE positive' = '#006D2C', 'TSE neutral' = '#f5d745', 'TSE negative' = '#FF0000')
color_mUCsubtype = c("Luminal_a" = "#33A02C", "Luminal_b" = "#1F78B4", "Stroma_rich" = "#ad3232", "Basal_squamous" = "#ff9900", "Non_specified" = "grey90")
color_consensusMIBC <- c("Basal/Squamous" = "#e0321b", "Stroma-rich" = "#f0d748",
                         "Luminal papillary" = "#7df0a7", "Luminal non-specified" = "#32a852",
                         "Luminal unstable" = "#327dfc", "Neuroendocrine-like" = "#f545d7")
color_targetedTherapy <- c("NECTIN4_amp" = '#d92b3c', "FGFR_mut" = "#1e73d4", "TSE_pos" = "#2ca324", "Other" = "grey80")


# List of phenotypic markers
genesPhenotypicMarkers <- c('CD44', 'CDH3', 'KRT1', 'KRT14', 'KRT16', 'KRT5', 'KRT6A', 'KRT6B', 'KRT6C', #Basal markers
                            'DSC1', 'DSC2', 'DSC3', 'DSG1', 'DSG2', 'DSG3', 'S100A7', 'S100A8', #squamous-differentiation markers
                            "COL1A1", "COL1A2", "COL6A1", "COL6A2", "COL6A3", "DCN", "THY1", "FAP",  # Stroma cell markers (CAF)
                            'CYP2J2', 'ERBB2', 'ERBB3', 'FGFR3', 'FOXA1', 'GATA3', 'GPX2', 'KRT18', 'KRT19', 'KRT20', 'KRT7', 'KRT8', 'PPARG', 'XBP1', 'UPK1A', 'UPK2', #Luminal
                            'CHGA', 'CHGB', 'SCG2', 'ENO2', 'SYP', 'NCAM1') #neuroendocrine markers
genesPhenotypicMarkers <- base::split(genesPhenotypicMarkers, c(rep("Basal", 9), rep("Squamous", 8), rep("Stroma", 8),
                                                                rep("Luminal", 16), rep("Neuroendocrine", 6)))


# Helper functions --------------------------------------------------------

annotationTheme <- function(){
  theme(legend.position = 'bottom', axis.ticks = element_blank(), axis.title.y = element_text(size = 8), axis.text.x = element_blank(), text=element_text(size=8, family='Helvetica'),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(NULL),
        plot.margin = ggplot2::margin(t = 0, 0,0,0, 'cm'),
        legend.margin = ggplot2::margin(t=-1, r=0, b=.5, l=0, unit='cm')
  )
}




# Get normalized RNA counts ----------------------------------------------------
rawCounts_DR314 <- round(read.csv(paste0(path.rna, "/rawCounts_full.txt"), sep="\t"))

# define colData
colData <- data.frame(row.names = colnames(rawCounts_DR314), sampleId = colnames(rawCounts_DR314))

# Normalize counts with DESeq2
DESeq2_DR314 <- DESeq2::DESeqDataSetFromMatrix(countData = rawCounts_DR314, colData = colData, design = ~1)
DESeq2_DR314_WT <- DESeq2::DESeq(DESeq2_DR314, parallel = T)

# complete normalization with vst
DESeq2_vst <- DESeq2::vst(DESeq2_DR314_WT, blind = T)

# get final normalized count Matrix
normalizedCountMatrix <- SummarizedExperiment::assay(DESeq2_vst)



# ================== calculate TSE score ==============
# Load the TSE classifier (available at: https://github.com/ANakauma/TSEscore_ICIs)
source("/TSEscore_ICIs/R/classifier/TSE_classify.R")
load("/TSEscore_ICIs/R/data/centroids_TSE.RData")

# apply function to the rna normalized count Matrix (rows = gene symbols, columns sampleId)
TSEclass <- TSE_classify(x = normalizedCountMatrix, centroids_TSE = centroids_TSE)
TSEclass <- TSEclass %>% dplyr::mutate(TSE_category = gsub("_", " ", TSE_category))


# ================== estimate mUC subtype ==============
# Load the mUC classififer (available at: https://bitbucket.org/ccbc/dr31_hmf_muc/)
source("/dr31_hmf_muc/classifier/mUC_classify.R")
load(file = "/dr31_hmf_muc/classifier/centroids_mUC.RData")

# apply function to the rna normalized count Matrix (rows = gene symbols, columns sampleId)
mUCsubtype <- classifymUC(x = normalizedCountMatrix, centroids_mUC = centroids_mUC)



# ================== estimate MIBC subtype ==============
patientMIBCclass <- consensusMIBC::getConsensusClass(data.frame(normalizedCountMatrix), gene_id = "hgnc_symbol")
patientMIBCclass$sampleId <- rownames(patientMIBCclass)

# give complete names to RNA subtypes
patientMIBCclass <- patientMIBCclass %>% dplyr::mutate(consensusSubtypeMIBC = ifelse(consensusClass == "Ba/Sq", "Basal/Squamous",
                                                                                     ifelse(consensusClass == "LumNS", "Luminal non-specified",
                                                                                            ifelse(consensusClass == "LumP", "Luminal papillary",
                                                                                                   ifelse(consensusClass == "LumU", "Luminal unstable",
                                                                                                          ifelse(consensusClass == "NE-like", "Neuroendocrine-like", "Stroma-rich"))))))





# ================== get NECTIN4 amplification status ==============
# read GISTIC CN (Gene names are outdated, be careful)
gisticFolder <- "/studycohorts/DR314/postHMF/GISTIC/"
gisticAllCNperGene_score <- readr::read_delim(paste0(gisticFolder, pattern = "all_data_by_genes.txt"),  delim = "\t") %>%
  reshape2::melt(id.vars = c("Gene Symbol", "Gene ID", "Cytoband")) %>%
  dplyr::select(geneName = `Gene Symbol`, Cytoband, sampleId = variable, gisticScore = value)
gisticAllCNperGene_cat <- readr::read_delim(paste0(gisticFolder, pattern = "all_thresholded.by_genes.txt"),  delim = "\t") %>%
  reshape2::melt(id.vars = c("Gene Symbol", "Locus ID", "Cytoband")) %>%
  dplyr::select(geneName = `Gene Symbol`, Cytoband, sampleId = variable, gisticCN = value)
gisticCNperGene <- gisticAllCNperGene_score %>% dplyr::left_join(gisticAllCNperGene_cat, by = c("geneName", "Cytoband", "sampleId"))

# select TACSTD2 and NECTIN4 (old name is PVRL4)
selGene_amp <- gisticCNperGene %>% dplyr::filter(geneName %in% c("PVRL4", "TACSTD2"))
selGene_amp <- selGene_amp %>% 
  dplyr::mutate(isGene_amp = ifelse(gisticCN == 2, "Yes", "No")) %>%
  dplyr::mutate(geneName = ifelse(geneName == "PVRL4", "NECTIN4", geneName)) %>%
  dplyr::mutate(geneName = factor(geneName, levels = c("TACSTD2", "NECTIN4")))



# ================== get mutations and fusions ==============

# Select mutations
selGene_mutations <- dplyr::filter(results.Cohort$combinedReport, SYMBOL %in% c("FGFR3", "FGFR2", "TACSTD2", "NECTIN4")) %>%
  dplyr::mutate(is_SNV_MNV_Indel_SV = ifelse(!is.na(Consequence.Mut) | !is.na(Consequence.SV), "Yes", "No")) %>%
  dplyr::mutate(SYMBOL = factor(SYMBOL, levels = c("NECTIN4", "TACSTD2", "FGFR2", "FGFR3")))

# get fusions
pathToFusions_list <- list.files(path=paste0(path.hmf, "HMF/somatic"), pattern = ".linx.fusion.tsv",
                                      recursive = TRUE, full.names = TRUE)
pathToFusions_list <- data.frame(pathToFusions_list = pathToFusions_list, sampleId = gsub("/linx.*", "", gsub(".*somatic/", "", pathToFusions_list)))
pathToFusions_list <- pathToFusions_list %>% dplyr::filter(sampleId %in% DR314.MetaData$sampleId)
  
# read tsv files and create a a table with all fusions
fusions_data_all <- pbapply::pblapply(pathToFusions_list$pathToFusions_list, function(iFile){
  fusionData <- read.table(file = iFile, sep = '\t', header = TRUE)
  if(nrow(fusionData) == 0){return(NULL)}
  fusionData$sampleId <- gsub("/linx.*", "", gsub(".*somatic/", "", iFile))
  return(fusionData)
}, cl = 5)

# make a table
fusions_data_all <- do.call(rbind.data.frame, fusions_data_all)

# select specific fusions (manually check gene names in table)
selGene_fusions <- fusions_data_all %>% dplyr::filter(grepl("FGFR3|FGFR2|NECTIN4", name)) %>%
  dplyr::mutate(mainGene = ifelse(grepl("FGFR3", name), "FGFR3", NA),
                mainGene = ifelse(grepl("FGFR2", name), "FGFR2", mainGene),
                mainGene = ifelse(grepl("NECTIN4", name), "NECTIN4", mainGene),
                isFused = "Yes") %>%
  dplyr::arrange(likelihood) %>%
  dplyr::distinct(mainGene, sampleId, .keep_all = TRUE) %>%
  dplyr::mutate(mainGene = factor(mainGene, levels = c("NECTIN4", "FGFR2", "FGFR3")))




###################################################################################################
################################################# Figure 1 ########################################
###################################################################################################


# ==================  Get all necessary information for biomarker stratification  ==============
matchDataToPlot <- dplyr::select(mUCsubtype, sampleId, mUC_subtype) %>%
  dplyr::left_join(dplyr::select(patientMIBCclass, sampleId, consensusSubtypeMIBC)) %>%
  dplyr::left_join(dplyr::select(TSEclass, sampleId, TSE_category)) %>%
  dplyr::left_join(dplyr::filter(selGene_amp, geneName == 'NECTIN4') %>% dplyr::select(sampleId, isNECTIN4_amp = isGene_amp)) %>%
  dplyr::left_join(dplyr::filter(selGene_mutations, SYMBOL == "FGFR3") %>% dplyr::select(sampleId = sample, FGFR3_mut = is_SNV_MNV_Indel_SV)) %>%
  dplyr::left_join(dplyr::filter(selGene_mutations, SYMBOL == "FGFR2") %>% dplyr::select(sampleId = sample, FGFR2_mut = is_SNV_MNV_Indel_SV)) %>%
  dplyr::left_join(dplyr::filter(selGene_fusions, mainGene == "FGFR3") %>% dplyr::select(sampleId, FGFR3_fusion = isFused)) %>%
  dplyr::left_join(dplyr::filter(selGene_fusions, mainGene == "FGFR2") %>% dplyr::select(sampleId, FGFR2_fusion = isFused)) %>%
  dplyr::left_join(dplyr::filter(selGene_fusions, mainGene == "NECTIN4") %>% dplyr::select(sampleId, NECTIN4_fusion = isFused)) %>%
  dplyr::mutate(FGFR3_mut_fusion = ifelse(FGFR3_mut == "Yes" | FGFR3_fusion == "Yes", 1, 0)) %>%
  dplyr::mutate(isFGFR_mut = ifelse(FGFR3_mut == "Yes" | FGFR3_fusion == "Yes" | FGFR2_mut == "Yes" | FGFR2_fusion == "Yes", 1, 0),
                isTSEpos = ifelse(TSE_category == "TSE positive", 1, 0)) %>%
  dplyr::mutate(isFGFR_mut = ifelse(is.na(isFGFR_mut), 0, isFGFR_mut),
                isTSEpos = ifelse(is.na(isTSEpos), 0, isTSEpos)) %>% 
  dplyr::mutate(groupTargetTherapy = ifelse(isNECTIN4_amp == "Yes", "NECTIN4_amp",
                                            ifelse(isFGFR_mut, "FGFR_mut",
                                                   ifelse(isTSEpos, "TSE_pos", "Other")))) %>%
  dplyr::mutate(groupTargetTherapy = factor(groupTargetTherapy, levels = c('NECTIN4_amp', "FGFR_mut", "TSE_pos", "Other")))



# Sort on mutually exclusiveness
memoData <- matchDataToPlot %>% dplyr::select(sampleId, isNECTIN4_amp, FGFR3_mut, FGFR2_mut, FGFR3_fusion, FGFR3_mut_fusion, mUC_subtype, TSE_category) %>%
  dplyr::mutate(isNECTIN4_amp = ifelse(isNECTIN4_amp == "No", 0, 1),
                FGFR3_mut = ifelse(is.na(FGFR3_mut), 0, 1),
                FGFR2_mut = ifelse(is.na(FGFR2_mut), 0, 1),
                FGFR3_fusion = ifelse(is.na(FGFR3_fusion), 0, 1),
                FGFR3_mut_fusion = ifelse(is.na(FGFR3_mut_fusion), 0, 1)) %>%
  dplyr::mutate(orderSamplesNum = ifelse(isNECTIN4_amp == 1, 10,
                                         ifelse(FGFR3_mut == 1, 9,
                                                ifelse(FGFR2_mut == 1, 8,
                                                       ifelse(FGFR3_fusion == 1, 7, 0))))) %>%
  dplyr::arrange(factor(mUC_subtype, levels = names(color_mUCsubtype))) %>%
  dplyr::arrange(factor(TSE_category, levels = names(color_TSEscore))) %>%
  dplyr::arrange(FGFR3_fusion) %>%
  dplyr::arrange(FGFR3_mut) %>%
  dplyr::arrange(FGFR2_mut) %>%
  dplyr::arrange(isNECTIN4_amp) %>%
  dplyr::arrange(-orderSamplesNum)
  

# Order axis
matchDataToPlot <- matchDataToPlot %>%
  dplyr::arrange(factor(sampleId, levels = memoData$sampleId))

# order samples from meta data
matchDataToPlot$sampleId <- factor(matchDataToPlot$sampleId, levels = matchDataToPlot$sampleId) 
matchDataToPlot$mUC_subtype <- factor(matchDataToPlot$mUC_subtype, levels = names(color_mUCsubtype))
selGene_amp$sampleId <- factor(selGene_amp$sampleId, levels = matchDataToPlot$sampleId)
selGene_fusions$sampleId <-  factor(selGene_fusions$sampleId, levels = matchDataToPlot$sampleId)
selGene_mutations$sample <- factor(selGene_mutations$sample, levels = matchDataToPlot$sampleId)



# plot NECTIN4 ampl (excl not matched RNA-DNA samples)
plot.Gene_ampl <- ggplot(dplyr::filter(selGene_amp, !is.na(sampleId)), aes(sampleId, y = geneName, fill = isGene_amp)) +
  geom_tile() +
  labs(y = "Gene\namplification", x = NULL) +
  annotationTheme() +
  scale_fill_manual('',
                    values = c('Yes' = '#02669c', 'No' = 'white'),
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))

# plot Fusions (excl not matched RNA-DNA samples)
plot.Gene_fusion <- ggplot(dplyr::filter(selGene_fusions, !is.na(sampleId)), aes(sampleId, y = mainGene, fill = isFused)) +
  geom_tile(na.rm = T) +
  labs(y = "Fusion\ngenes", x = NULL) +
  scale_x_discrete(drop=FALSE) +
  annotationTheme() +
  scale_fill_manual('',
                    values = c('Yes' = '#02669c', 'No' = 'white'), 
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))



# plot NECTIN4 ampl (excl not matched RNA-DNA samples)
plot.Gene_nonSynMutations <- ggplot(dplyr::filter(selGene_mutations, !is.na(sample)), aes(sample, y = SYMBOL, fill = is_SNV_MNV_Indel_SV)) +
  geom_tile(na.rm = F) +
  labs(y = "Non-syn\nmutations", x = NULL) +
  scale_x_discrete(drop=FALSE) +
  annotationTheme() +
  scale_fill_manual('',
                    values = c('Yes' = '#02669c', 'No' = 'white'), 
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))



# plot TSE score
plot.TSEscore <- ggplot(matchDataToPlot, aes(sampleId, y = "TSE score", fill = TSE_category)) +
  geom_tile() +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('TSE score',
                    values = color_TSEscore,
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot.TSEscore <- cowplot::get_legend(
  plot.TSEscore + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.TSEscore <- plot.TSEscore + theme(legend.position = 'none')



# plot mUC subtypes
plot.mUCsubtype <- ggplot(matchDataToPlot, aes(sampleId, y = "mUC subtype", fill = mUC_subtype)) +
  geom_tile() +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('mUC subtype',
                    values = color_mUCsubtype,
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot.mUCsubtype <- cowplot::get_legend(
  plot.mUCsubtype + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.mUCsubtype <- plot.mUCsubtype + theme(legend.position = 'none')



# plot consensus MIBC subtypes
plot.MIBCsubtype <- ggplot(matchDataToPlot, aes(sampleId, y = "Consensus MIBC subtype", fill = consensusSubtypeMIBC)) +
  geom_tile() +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('MIBC subtype',
                    values = color_consensusMIBC,
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'right',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot.MIBCsubtype <- cowplot::get_legend(
  plot.MIBCsubtype + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot.MIBCsubtype <- plot.MIBCsubtype + theme(legend.position = 'none')









############################### Plot expression of selected genes
selectedGenes_exp <- c("NECTIN4", "FGFR3", "FGFR2", "FGFR1", 'TACSTD2', "GATA3", "PPARG", "CTLA4", "CD274", "PDCD1")
# Get only selected samples and genes
normMatrixRNA_selGenes <- normalizedCountMatrix[, colnames(normalizedCountMatrix) %in% matchDataToPlot$sampleId]
normMatrixRNA_selGenes <- normMatrixRNA_selGenes[rownames(normMatrixRNA_selGenes) %in% selectedGenes_exp,]
normMatrixRNA_selGenes  <- normMatrixRNA_selGenes - rowMedians(normMatrixRNA_selGenes)

# Get data frame of expression matrix
normMatrixRNA_selGenes <- as.data.frame(t(normMatrixRNA_selGenes))
normMatrixRNA_selGenes$sampleId <- rownames(normMatrixRNA_selGenes)

# melt data
normMatrixRNA_selGenes <- reshape2::melt(normMatrixRNA_selGenes, id.vars = c("sampleId"))

# order data
normMatrixRNA_selGenes <- normMatrixRNA_selGenes %>%
  dplyr::mutate(variable = factor(variable, levels = rev(selectedGenes_exp)),
                sampleId = factor(sampleId, levels = matchDataToPlot$sampleId))

# change max min expression value
minExpValue <- 2.2
normMatrixRNA_selGenes <- normMatrixRNA_selGenes %>%
  dplyr::mutate(value_tmp = ifelse(value < -minExpValue, -minExpValue,
                                   ifelse(value > minExpValue, minExpValue, value)))

plot_selGeneExpression <- ggplot(normMatrixRNA_selGenes, aes(x = sampleId, y = variable, fill = value_tmp)) +
  geom_tile(size = 0.25, na.rm = F, alpha = 1, height = 1, width = 1) +
  scale_fill_gradient2(low = "green", mid = 'black', high = 'red',
                       breaks=c(-2, 0, 2),labels=c("<-2", 0, ">2"),
                       guide = guide_colorbar(title = 'Normalized expression',
                                              title.position = 'top', title.hjust = 0.5, barwidth = 5, barheight = 0.5)) +
  annotationTheme() +
  labs(x = NULL, y = "Gene\nexpression") +
  theme(
    legend.position = 'bottom',
    text=element_text(size=8, family='Helvetica'),
    axis.text=element_text(size=8),
    axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))

# extract legends
legend_plot_selGeneExpression <- cowplot::get_legend(
  plot_selGeneExpression + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_selGeneExpression <- plot_selGeneExpression + theme(legend.position = 'none')






# Plot phenotypic signatures ------------------------------------------------
gsvapar <- GSVA::gsvaParam(normalizedCountMatrix, genesPhenotypicMarkers)
results_phenotypicMarkers <- GSVA::gsva(gsvapar)

# Melt results
results_phenotypicMarkers <- as.data.frame(results_phenotypicMarkers)
results_phenotypicMarkers$phenotype <- rownames(results_phenotypicMarkers)
results_phenotypicMarkers <- reshape2::melt(results_phenotypicMarkers, id.vars = "phenotype")

# only keep selected samples
results_phenotypicMarkers <- results_phenotypicMarkers %>% dplyr::filter(variable %in% matchDataToPlot$sampleId)

# order samples and phenotypes
results_phenotypicMarkers <- results_phenotypicMarkers %>%
  dplyr::mutate(variable = factor(variable, levels = matchDataToPlot$sampleId),
                phenotype = factor(phenotype, levels = rev(c("Luminal", "Stroma" , "Basal", "Squamous", "Neuroendocrine"))))


# Set max.min values
minExp <- 0.55
results_phenotypicMarkers  <- results_phenotypicMarkers %>%
  dplyr::mutate(value_tmp = ifelse(value > minExp, minExp, ifelse(value < -minExp, -minExp, value)))

# plot phenotypic score
plot_phenotypicMarkers <- ggplot(results_phenotypicMarkers, aes(x = variable, y = phenotype, fill = value_tmp)) +
  geom_tile(size = 0.25, na.rm = F, alpha = 1, height = 1, width = 1) +
  scale_fill_gradient2(low = "#3288BD", mid = 'white', high = '#9E0142',
                       breaks=c(-0.5, 0, 0.5),labels=c("<-0.5", 0, ">0.5"),
                       guide = guide_colorbar(title = 'Phenotypic signature',
                                              title.position = 'top', title.hjust = 0.5, barwidth = 5, barheight = .5)) +
  annotationTheme() +
  labs(x = NULL, y = 'Phenotypic\nsignature') +
  theme(
    legend.position = 'bottom',
    text=element_text(size=8, family='Helvetica'),
    axis.text=element_text(size=8),
    axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'))


# extract legends
legend_plot_phenotypicMarkers <- cowplot::get_legend(
  plot_phenotypicMarkers + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# delete legend from main plot
plot_phenotypicMarkers <- plot_phenotypicMarkers + theme(legend.position = 'none')







# Plot genomic Landscape Figure 1A ---------------------------------------------------

# Plot sample-level overview.
samplePlot <- cowplot::plot_grid(plot.Gene_ampl,
                                 plot.Gene_nonSynMutations,
                                 plot.Gene_fusion,
                                 plot.TSEscore,
                                 plot.mUCsubtype,
                                 plot.MIBCsubtype,
                                 plot_selGeneExpression,
                                 plot_phenotypicMarkers,
                                 ncol = 1, align = 'v', axis = 'tblr',
                                 rel_heights = c(0.5, 0.7, 0.7, 0.3, 0.3, 0.3, 2.1, 1.3)
)

legend_samplePlot <- cowplot::plot_grid(legend_plot.TSEscore,
                                        legend_plot.mUCsubtype,
                                        legend_plot.MIBCsubtype,
                                        legend_plot_selGeneExpression,
                                        legend_plot_phenotypicMarkers,
                                        nrow = 1, align = 'v')


pdf(paste0(odir,"/subtype_overview.pdf"), width = 8, height = 4.3)#, width = 14, height = 21)
cowplot::plot_grid(samplePlot, legend_samplePlot, ncol = 1, rel_heights = c(1, 0.2), align = 'h', axis = 'tblr')
dev.off()






# Mutually exclusive test ========================================

# Sort on mutually exclusiveness.
matchDataToPlot_summary <- matchDataToPlot %>% dplyr::select(sampleId, isNECTIN4_amp, isFGFR_mut,  isTSEpos) %>%
  dplyr::mutate(isNECTIN4_amp = ifelse(isNECTIN4_amp == "Yes", 1, 0))
rownames(matchDataToPlot_summary) <- matchDataToPlot_summary$sampleId
matchDataToPlot_summary$sampleId <- NULL
matchDataToPlot_summary <- t(matchDataToPlot_summary)

# Use Rediscover R package to find exclusivity p-values
PMA <- Rediscover::getPM(as.matrix(matchDataToPlot_summary))
exclusivity_test <- as.matrix(Rediscover::getMutex(A=as.matrix(matchDataToPlot_summary), PM=PMA, lower.tail = TRUE))
colnames(exclusivity_test) <- rownames(matchDataToPlot_summary)
rownames(exclusivity_test) <- rownames(matchDataToPlot_summary)
exclusivity_test <- -log10(exclusivity_test)
exclusivity_test[row(exclusivity_test)==col(exclusivity_test)] <- 0

# plot p-values of mutual exclusivity
plot_mutalExcl <- function() {
  corrplot::corrplot(exclusivity_test, type = "full", diag = FALSE, is.corr = FALSE, tl.col = "black",
                     col = colorRampPalette(c("white", "#f78f8f", "red", "#6D0000"))(8), col.lim = c(0, 4), tl.cex = 0.6, tl.pos = 'lt')
}


# export results
pdf(paste0(odir, "Fig1b_excTest.pdf"), width = 7, height = 3.5)
plot_mutalExcl()
dev.off()





############################### Compare gene expression between patient groups
selectedGenes_exp_stats <- c("NECTIN4", "FGFR3", "FGFR2", "CTLA4", "CD274", "PDCD1")

# Get only selected samples and genes
normMatrixRNA_selGenes <- normalizedCountMatrix[, colnames(normalizedCountMatrix) %in% matchDataToPlot$sampleId]
normMatrixRNA_selGenes <- normMatrixRNA_selGenes[rownames(normMatrixRNA_selGenes) %in% selectedGenes_exp_stats,]
normMatrixRNA_selGenes  <- normMatrixRNA_selGenes - rowMedians(normMatrixRNA_selGenes)

# Get data frame of expression matrix
normMatrixRNA_selGenes <- as.data.frame(t(normMatrixRNA_selGenes))
normMatrixRNA_selGenes$sampleId <- rownames(normMatrixRNA_selGenes)

# melt data
normMatrixRNA_selGenes <- reshape2::melt(normMatrixRNA_selGenes, id.vars = c("sampleId"))

# Define group with no NECTIN4 amp, FGFR2/3 mut, TSE positive and others
normMatrixRNA_selGenes_Other <-  normMatrixRNA_selGenes %>% dplyr::left_join(dplyr::select(matchDataToPlot, sampleId, isNECTIN4_amp, isFGFR_mut, isTSEpos)) %>%
  dplyr::mutate(patientGroup = ifelse(isNECTIN4_amp == "No" & isFGFR_mut == 0 & isTSEpos == 0, "Other", "Mut")) %>%
  dplyr::filter(patientGroup == "Other") %>%
  dplyr::mutate(isNECTIN4_amp = NULL, isFGFR_mut = NULL, isTSEpos = NULL)

# add info on mutations and TSE positivity
normMatrixRNA_selGenes_nectin4 <- normMatrixRNA_selGenes %>% dplyr::left_join(dplyr::select(matchDataToPlot, sampleId, patientGroup = isNECTIN4_amp)) %>%
  dplyr::mutate(patientGroup = ifelse(patientGroup == "Yes", "NECTIN4 amp", "Other"))
normMatrixRNA_selGenes_fgfr <- normMatrixRNA_selGenes %>% dplyr::left_join(dplyr::select(matchDataToPlot, sampleId, patientGroup = isFGFR_mut)) %>%
  dplyr::mutate(patientGroup = ifelse(patientGroup == 1, "FGFR2/3 mut", "Other"))
normMatrixRNA_selGenes_TSEpos <- normMatrixRNA_selGenes %>% dplyr::left_join(dplyr::select(matchDataToPlot, sampleId, patientGroup = isTSEpos)) %>%
  dplyr::mutate(patientGroup = ifelse(patientGroup == 1, "TSE positive", "Other"))
normMatrixRNA_selGenes <- rbind(normMatrixRNA_selGenes_nectin4, normMatrixRNA_selGenes_fgfr, normMatrixRNA_selGenes_TSEpos) %>%
  dplyr::filter(patientGroup %in% c("NECTIN4 amp", "FGFR2/3 mut", "TSE positive"))

# add "Other"
normMatrixRNA_selGenes <- rbind(normMatrixRNA_selGenes, normMatrixRNA_selGenes_Other)


# order data
normMatrixRNA_selGenes <- normMatrixRNA_selGenes %>%
  dplyr::mutate(variable = factor(variable, levels = c("NECTIN4", "FGFR3", "FGFR2", "CTLA4", "CD274", "PDCD1")),
                patientGroup = factor(patientGroup, levels = c("NECTIN4 amp", "FGFR2/3 mut", "TSE positive", "Other")))


# stats
#  Perform kruskal test groups
geneExp_stat <- normMatrixRNA_selGenes %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(p.value = stats::kruskal.test(value ~ patientGroup)$p.value) %>% dplyr::ungroup()
geneExp_stat$pAdj <- p.adjust(geneExp_stat$p.value, method = "BH")
geneExp_stat <- geneExp_stat %>%
  dplyr::mutate(pAdj.label = ifelse(pAdj < 0.001, 'q<0.001',
                                    ifelse(pAdj < 0.01, paste0("q=", round(pAdj, 3)), paste0("q=", round(pAdj, 2)))))

# Plot expression genes between groups
plot_selGeneExperssion_stats <- ggplot(normMatrixRNA_selGenes, aes(x = variable, y = value, fill = patientGroup)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", linewidth = 0.4) +
  geom_boxplot(notch = F, alpha = .8, size = 0.3, width=0.85, outlier.shape = NA) +
  geom_point(aes(x = variable, fill = patientGroup),  alpha = 0.2, cex = 0.5, stroke = 0.2,
             shape = 21, position = position_jitterdodge(jitter.width = 0.15)) +
  ylab('Gene expression\n(normalized counts)') +
  scale_fill_manual(na.translate=FALSE,
                    guide = guide_legend(title = NULL, title.position = 'top',
                                         title.hjust = 0.6, nrow = 1, keywidth = 1, keyheight = 1),
                    name = NULL,
                    values = c('NECTIN4 amp' = '#d92b3c', 'FGFR2/3 mut' = '#1e73d4', 'TSE positive' = '#2ca324', 'Other' = 'grey80')) +
  # Add adjusted p-values from kruskal.test
  geom_text(data=geneExp_stat, aes(x = c(1, 2, 3, 4, 5, 6), y = rep(6, 6),
                                   label=pAdj.label), size = 2.5, inherit.aes = FALSE) +
  theme(
    legend.position = 'bottom',
    axis.title.y = element_text(size = 8),
    text=element_text(size=8, family='Helvetica'),
    panel.spacing=unit(0, "lines"),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x = element_text(size = 8, angle=30, hjust=1, vjust=1.0),
    axis.text.y = element_text(size = 8),
    panel.background = element_rect(fill = 'white', colour = NA),
    panel.border = element_rect(fill = NA, colour = 'grey20'),
    strip.background = element_rect(colour = 'grey20', fill = 'white')
  )


# Export results
pdf(paste0(odir, "FigS1c_stats_HMF.pdf"), width = 3.25, height = 2)
plot_selGeneExperssion_stats
dev.off()





# Plot expression genes and pair-wise comparison between groups
plot_selGeneExperssion_stats_extra <- ggplot(normMatrixRNA_selGenes, aes(x = patientGroup, y = value, fill = patientGroup)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey", linewidth = 0.4) +
  geom_boxplot(notch = F, alpha = .8, size = 0.3, width=0.85, outlier.shape = NA) +
  geom_point(aes(x = patientGroup, fill = patientGroup),  alpha = 0.2, cex = 0.5, stroke = 0.2,
             shape = 21, position = position_jitterdodge(jitter.width = 0.15)) +
  facet_wrap(~variable, nrow = 1) +
  ylab('Gene expression\n(normalized counts)') +
  scale_fill_manual(na.translate=FALSE,
                    guide = guide_legend(title = NULL, title.position = 'top',
                                         title.hjust = 0.6, nrow = 1, keywidth = 1, keyheight = 1),
                    name = NULL,
                    values = c('NECTIN4 amp' = '#d92b3c', 'FGFR2/3 mut' = '#1e73d4', 'TSE positive' = '#2ca324', 'Other' = 'grey80')) +
  theme(
  legend.position = 'bottom',
  axis.title.y = element_text(size = 8),
  text=element_text(size=8, family='Helvetica'),
  panel.spacing=unit(0, "lines"),
  axis.title.x=element_blank(),
  axis.ticks.x=element_blank(),
  axis.text.x = element_text(size = 8, angle=30, hjust=1, vjust=1.0),
  axis.text.y = element_text(size = 8),
  panel.background = element_rect(fill = 'white', colour = NA),
  panel.border = element_rect(fill = NA, colour = 'grey20'),
  strip.background = element_rect(colour = 'grey20', fill = 'white')
) +
  ggsignif::stat_signif(show.legend = F,
                        comparisons=list(c('NECTIN4 amp', 'FGFR2/3 mut'), c('TSE positive', 'Other'),
                                         c('FGFR2/3 mut', 'TSE positive'),
                                         c('NECTIN4 amp', 'TSE positive'), c('NECTIN4 amp', 'Other'),
                                         c('FGFR2/3 mut', 'Other')),
                        test = 'wilcox.test',
                        test.args=list(paired=FALSE),
                        map_signif_level = F,
                        step_increase = .04,
                        textsize = 2.5,
                        color = 'black', tip_length = .01)

# Export results
pdf(paste0(odir, "FigS1c_stats_HMF_tmp.pdf"), width = 3.25, height = 2.5)
plot_selGeneExperssion_stats_extra
dev.off()








###################################################################################################
################################################# Figure 2 ########################################
###################################################################################################



# Plot summary of treatment options and markers overlap ---------------------------------------------------

# Sort data according to NECTIN4 amp, FGFR2/3 muts and TSE+
matchDataToPlot <- matchDataToPlot %>%
  dplyr::arrange(-isFGFR_mut, -isTSEpos) %>% dplyr::arrange(factor(isNECTIN4_amp, levels = c("Yes", "No")))
matchDataToPlot$sampleId <- factor(matchDataToPlot$sampleId, levels = matchDataToPlot$sampleId) 

# plot targeted therapy groups
plot.targetedTherapy <- ggplot(matchDataToPlot, aes(sampleId, y = "Patient group", fill = groupTargetTherapy)) +
  geom_tile() +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('',
                    values = color_targetedTherapy,
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))



# plot targeted therapy groups
plot.NECTIN4amp <- ggplot(matchDataToPlot, aes(sampleId, y = "NECTIN4 amp", fill = isNECTIN4_amp)) +
  geom_tile() +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('',
                    values = c('Yes' = '#d92b3c', 'No' = 'white'), 
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))


# plot targeted therapy groups
plot.FGFRmuts <- ggplot(matchDataToPlot, aes(sampleId, y = "FGFR2/3 mut", fill = as.factor(isFGFR_mut))) +
  geom_tile() +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('',
                    values = c('1' = "#1e73d4", '0' = 'white'), 
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))




# plot targeted therapy groups
plot.TSEpos <- ggplot(matchDataToPlot, aes(sampleId, y = "+TSE score", fill = as.factor(isTSEpos))) +
  geom_tile() +
  labs(y = NULL, x = NULL) +
  annotationTheme() +
  scale_fill_manual('',
                    values = c('1' = "#2ca324", '0' = 'white'), 
                    guide = guide_legend(title.position = 'top', title.hjust = 0, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.y = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 8),
        text=element_text(size=8, family='Helvetica', colour = "black"),
        panel.border = element_rect(fill = NA, colour = 'grey20'))






#--------------------------- plot Survival -------------------------------------------
# add treatment groups 
survival_DR314 <- DR314.MetaData %>%
  dplyr::left_join(dplyr::select(matchDataToPlot, sampleId, groupTargetTherapy)) %>%
  dplyr::filter(!is.na(groupTargetTherapy)) %>%
  dplyr::mutate(biopsySite = ifelse(biopsySite %in% c('Liver', 'Lymph node'), biopsySite, "Other"),
                primaryTumor = ifelse(primaryTumorSubLocation %in% c('Bladder', 'UTUC'), primaryTumorSubLocation, 'Other'))
  

# Cox regression analysis -------------------------------

# plot ICI subgroup 
coxRegressionAnalysis <-  dplyr::filter(survival_DR314, receivedImmuno == 1) %>%
  survivalAnalysis::analyse_multivariate(vars(OS_months_treatmentDate, isDead),
                                         covariates = vars(groupTargetTherapy,
                                                           biopsySite,
                                                           primaryTumor,
                                                           treatmentGroup,
                                                           ageAtBiopsy,
                                                           gender,),
                                         reference_level_dict=c(gender="male", biopsySite='Other', primaryTumor='Other'))

# save results
pdf(paste0(odir, "MVA_DR314_immuno.pdf"), width = 5, height = 2.5)
survivalAnalysis::forest_plot(coxRegressionAnalysis,
                              endpoint_labeller = c(OS_months_treatmentDate="OS"),
                              labels_displayed = c("endpoint", "factor", "n"),
                              ggtheme = ggplot2::theme_bw(base_size = 8),
                              relative_widths = c(1, 1.5, 1),
                              HR_x_limits = c(0.1, 16.5),
                              HR_x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8))
dev.off()



# plot chemo subgroup
coxRegressionAnalysis <- dplyr::filter(survival_DR314, receivedChemo == 1) %>%
  survivalAnalysis::analyse_multivariate(vars(OS_months_treatmentDate, isDead),
                                         covariates = vars(groupTargetTherapy,
                                                           biopsySite,
                                                           primaryTumor,
                                                           ageAtBiopsy,
                                                           gender),
                                         reference_level_dict=c(gender="male", biopsySite='Other', primaryTumor='Other'))

# save results
pdf(paste0(odir, "MVA_DR314_chemo.pdf"), width = 5, height = 2.5)
survivalAnalysis::forest_plot(coxRegressionAnalysis,
                              endpoint_labeller = c(OS_months_treatmentDate="OS"),
                              labels_displayed = c("endpoint", "factor", "n"),
                              ggtheme = ggplot2::theme_bw(base_size = 8),
                              relative_widths = c(1, 1.5, 1),
                              HR_x_limits = c(0.05, 16.5),
                              HR_x_breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8))
dev.off()




# Survival plots -------------------------------
# get survival model ******
fit_OS_immuno <- survminer::surv_fit(survival::Surv(OS_months_treatmentDate,  isDead) ~ groupTargetTherapy, data = dplyr::filter(survival_DR314, receivedImmuno == 1))
fit_OS_chemo <- survminer::surv_fit(survival::Surv(OS_months_treatmentDate,  isDead) ~ groupTargetTherapy, data = dplyr::filter(survival_DR314, receivedChemo == 1))

# survival plot for ICI treated patients
plot_OSCurve_Therapy_immuno <- survminer::ggsurvplot(fit_OS_immuno, data =  dplyr::filter(survival_DR314, receivedImmuno == 1),
                                                     size = 1,                 # change line size
                                                     palette =  c("grey80", "#2ca324", "#1e73d4", '#d92b3c'),
                                                     xlab = "Months",
                                                     ylab = "OS probability",
                                                     xlim = c(0,24),
                                                     legend.title = "Targeted therapy group",
                                                     legend = c(0.7,0.7),
                                                     pval.method.size = 5,
                                                     break.time.by = 6,
                                                     pval = TRUE,              # Add p-value (default log-rank)
                                                     pval.size = 3,
                                                     pval.coord = c(2.5, 0.1),
                                                     risk.table = TRUE,        # Add risk table
                                                     risk.table.fontsize = 2.5,
                                                     tables.y.text.col = FALSE, 
                                                     legend.labs = levels(survival_DR314$groupTargetTherapy),    # Change legend labels
                                                     ggtheme = theme(
                                                       legend.position = 'bottom',
                                                       text=element_text(size=10, family='Helvetica'),
                                                       plot.title = element_text(size=8, hjust = 0),
                                                       legend.title=element_text(size=8, hjust = 0.5), 
                                                       legend.text=element_text(size=8),
                                                       legend.key = element_rect(color = NA, fill = NA),
                                                       legend.key.size = unit(0.25, "cm"),
                                                       legend.spacing.x = unit(0.1, "line"),
                                                       legend.spacing.y = unit(0, "line"),
                                                       legend.background	= element_rect(fill = NA, colour = NA),
                                                       panel.grid.major.x = element_blank(),
                                                       panel.grid.minor.x = element_blank(),
                                                       panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(),
                                                       panel.background = element_rect(fill = 'white', colour = NA),
                                                       panel.border = element_rect(fill = NA, colour = 'grey20'))
)      # Change ggplot2 theme)



# survival plot for chemo treated patients
plot_OSCurve_Therapy_chemo <- survminer::ggsurvplot(fit_OS_chemo, data =  dplyr::filter(survival_DR314, receivedChemo == 1),
                                                    size = 1,                 # change line size
                                                    palette =  c("grey80", "#2ca324", "#1e73d4", '#d92b3c'),
                                                    xlab = "Months",
                                                    ylab = "OS probability",
                                                    xlim = c(0,24),
                                                    legend.title = "Targeted therapy group",
                                                    legend = c(0.7,0.7),
                                                    pval.method.size = 5,
                                                    break.time.by = 6,
                                                    pval = TRUE,              # Add p-value (default log-rank)
                                                    pval.size = 3,
                                                    pval.coord = c(2.5, 0.1),
                                                    risk.table = TRUE,        # Add risk table
                                                    risk.table.fontsize = 2.5,
                                                    tables.y.text.col = FALSE, 
                                                    legend.labs = levels(survival_DR314$groupTargetTherapy),    # Change legend labels
                                                    ggtheme = theme(
                                                      legend.position = 'bottom',
                                                      text=element_text(size=10, family='Helvetica'),
                                                      plot.title = element_text(size=8, hjust = 0),
                                                      legend.title=element_text(size=8, hjust = 0.5), 
                                                      legend.text=element_text(size=8),
                                                      legend.key = element_rect(color = NA, fill = NA),
                                                      legend.key.size = unit(0.25, "cm"),
                                                      legend.spacing.x = unit(0.1, "line"),
                                                      legend.spacing.y = unit(0, "line"),
                                                      legend.background	= element_rect(fill = NA, colour = NA),
                                                      panel.grid.major.x = element_blank(),
                                                      panel.grid.minor.x = element_blank(),
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),
                                                      panel.background = element_rect(fill = 'white', colour = NA),
                                                      panel.border = element_rect(fill = NA, colour = 'grey20'))
)      # Change ggplot2 theme)



# save figure OS for biomarkers-guided therapy stratification
pdf(paste0(odir,"Survival_HMF_immuno.pdf"),width = 4, height = 4.5)
cowplot::plot_grid(
  plot.targetedTherapy,
  plot.NECTIN4amp,
  plot.FGFRmuts,
  plot.TSEpos,
  plot_OSCurve_Therapy_immuno$plot,
  plot_OSCurve_Therapy_immuno$table,
  align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(0.1, 0.07, 0.07, 0.07, 1, 0.4))
dev.off()



# save figure OS for biomarkers-guided therapy stratification
pdf(paste0(odir,"Survival_HMF_chemo.pdf"),width = 4, height = 4.5)
cowplot::plot_grid(
  plot.targetedTherapy,
  plot.NECTIN4amp,
  plot.FGFRmuts,
  plot.TSEpos,
  plot_OSCurve_Therapy_chemo$plot,
  plot_OSCurve_Therapy_chemo$table,
  align = 'v', axis = 'tblr', ncol = 1, rel_heights = c(0.1, 0.07, 0.07, 0.07, 1, 0.4))
dev.off()




 