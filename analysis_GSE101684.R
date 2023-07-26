# CircRNA - GSE101684 NSCLC
# Arraystar Human circRNA array V2
# 4 normal samples + 4 tumor samples

## Packages
library(limma) #microarray analysis
library(Biobase) #EList to ExpressionSet conversion
library(arrayQualityMetrics) #QC
library(ggplot2) #PCA
library(ggrepel) #PCA
library(tidyr) #relative log expression
library(dendextend) #dendrograms
library(quantro) #normalization method selection
library(pvca) #batch effect visualization
library(plyr) #matrix matches search
library(xlsx) #data export
library(EnhancedVolcano) #volcano plots
library(gplots) #heatmap

## Custom functions ----
source("./customFunctions.R")

## general variables
my_colors <- c(rep("yellow", 4), rep("chartreuse3", 4))

## File list ----
#*make a file list for the analysis
files <- list.files(path="./RawData")
#*export file names to txt file for design matrix
write.table(files, file="samples.txt", col.names="FileName", row.names=FALSE)
#*check file list
write.table(files, col.names="FileName")

## Data import ----
#*expression data
dat <- read.maimages(files, "agilent.median", green.only=TRUE, path="./RawData")
#*check imported data (two numbers, first is no. of probes, second is no. of arrays)
dim(dat) 
#*phenoData
targets <- readTargets(file="samples_phenoData.txt")
targets
#*factors for each type of phenoData
Pair <- factor(targets$Pair)
Group <- factor(targets$Group, levels = c("Tumor", "Normal"))
Sex <- factor(targets$Sex, levels=c("Male", "Female"))

summary(Pair)
summary(Group)
summary(Sex)

#*phenoData for ExpressionSet data (needed for arrayQualityMetrics and PVCA)
annotated_targets <- AnnotatedDataFrame(targets)

## Annotation ----
#*annotation import
probeAnnotation <- read.columns("./Annotation/complete_annotations_V2.txt")
#*search for matches, attach circRNA names to dat variable
dat$genes$circRNA <- probeAnnotation[match(dat$genes$ProbeName, probeAnnotation$ProbeName),2]
head(dat) #check circRNA names

## Raw data quality control ----
#*expression data (E, Eb)
dat_E <- dat$E
dat_Eb <- dat$Eb

#*boxplot - check signal intensities
#//las - twisted lables, range - whisker length 
#//0-all data
boxplot(log2(dat_E), range=0, las=3, ylab="log2 intensity", main="Raw data", 
        col=my_colors)
boxplot(log2(dat_Eb), range=0, las=3, ylab="log2 intensity", main="Raw data background", 
        col=my_colors)

#*plotDensities
plotDensities(log2(dat$E), legend="topright", main="Raw data",
              Group=targets$Group)
plotDensities(log2(dat$Eb), legend="topright", main="Raw data background",
              Group=targets$Group)

#*PCA
plotCustomPCA(log2(dat_E), labels=targets$Name, groups=targets$Group,
              title="Raw data", colors=c("green", "yellow", "red", "blue"),
              size=4)

#*relative log expression
plotCustomRLE(log2(dat_E), colors=my_colors, title="raw", ylim = c(-3, 3))

#*dendrogram
plotCustomDendrogram(input_data = dat_E, title = "Raw data", groups = targets$Group, 
                     names = targets$Name, colors = c(Tumor="deeppink2", Normal="deepskyblue3"), 
                     groups_2 = targets$Sex, colors_2 = c(M="green3", F="darkorange1"), 
                     bar_label = "Sex", hanging=-0.08)
legend("topright", title="Legend", fill = c("deeppink2", "deepskyblue3", "green3", "darkorange1"),
       legend = as.character(c("Tumor", "Normal", "Male", "Female")))

#*ArrayQualityMetrics
arrayQualityMetrics(expressionset=ExpressionSet(assayData=dat_E, phenoData=annotated_targets), 
                    #data - dat_E, converted to ExpressionSet
                    outdir=file.path("./", "QC_Raw_E"), #folder destination
                    force=TRUE, #forces empty folder if already exists
                    do.logtransform=TRUE, #log2 data transformation (you need that until normalization)
                    intgroup="Group") #coloring according to group
arrayQualityMetrics(expressionset=ExpressionSet(assayData=dat_Eb, phenoData=annotated_targets),
                    outdir=file.path("./", "QC_Raw_Eb"), 
                    force=TRUE, do.logtransform=TRUE, intgroup="Group")
dev.off() #exit folder

## Background correction + quality control ----
#*Foreground - background plots 
#create or overwrite folder
dir.create("FB_plots", recursive = TRUE, showWarnings = FALSE)
for (i in c(1:length(files))){
        directory <- paste(c("FB_plots", "/"), collapse = "")
        file1 <-  paste(c(paste(c("Array", i), collapse = ""), "png"), collapse = ".")
        filename <- paste(directory, file1, collapse = "")
        png(file=filename)
        plotFB(dat, array = i)
        dev.off()
}

#*Background correction
#//methods: "auto", "none", "subtract", "half", "minimum", "movingmin", "edwards", "normexp"
#//
dat_bg <- backgroundCorrect(dat, method="movingmin", offset = 1)

#*plotDensities intensities of arrays (to check for differences between arrays)
#//legend="topright"
plotDensities(dat_bg, legend="topright", main="Background corrected data",
              Group=targets$Group)

#*boxplot - check signal intensities
#//las - twisted lables, range - whisker length 
#//0-all data
boxplot(log2(dat_bg$E), range=0, las=3, ylab="log2 intensity", main="Background corrected data",
        col=my_colors)

#*relative log expression
plotCustomRLE(log2(dat_bg$E), colors=my_colors, title="background corrected",
              ylim = c(-2, 3))

## Data normalization ----
#//methods (for normalizeBetweenArrays): "none", "scale", "quantile", "cyclicloess"
#//for two-color arrays: previous + "Aquantile", "Gquantile", "Rquantile", "Tquantile"

#*select normalization method (using quantro package)
#//very important if plotDensites shows differences between groups
#//if differences are due to biological variability, quantile normalization is not appropriate!
#//use cyclicloess instead
qtest <- quantro(object = dat_bg$E, groupFactor = Group)
qtest
summary(qtest)
anova(qtest)
quantroStat(qtest)
qtestPerm <- quantro(object = dat_bg$E, groupFactor = Group, B=1000)
qtestPerm
quantroPlot(qtestPerm)

#*normalization
dat_norm_0 <- normalizeBetweenArrays(dat_bg, method="cyclicloess")
#average duplicated probes (more common in V1 arrays)
dat_norm <- avereps(dat_norm_0, ID=dat_norm_0$genes$ProbeName)

## Normalized data quality control ----
#*expression data
dat_norm_E <- dat_norm$E

#*plotDensities with normalized data should match between arrays
plotDensities(dat_norm, legend="topright", main="Normalized data",
              Group=targets$Group)

#*boxplot - check signal intensities
#//las - twisted lables, range - whisker length 
#//0-all data
boxplot(dat_norm$E, range=0, las=3, ylab="log2 intensity", main="Normalized data", 
        col=my_colors)

#*PCA
plotCustomPCA(dat_norm_E, labels=targets$Name, groups=targets$Group,
              title="Normalized data", colors=c("green", "yellow", "red", "blue"),
              size=4)

#*batch effect detection
#graph
plotCustomPVCA(input_data = dat_norm$E, pheno_data = annotated_targets, 
               batch_factors = c("Group", "Sex"))

#*dendrogram
plotCustomDendrogram(input_data = dat_norm_E, title = "Normalized data", names = targets$Name,
                     groups = targets$Group, colors = c(Tumor="deeppink2", Normal="deepskyblue3"), 
                     groups_2 = targets$Sex, colors_2 = c(M="green3", F="darkorange1"),
                     bar_label = "Sex", hanging = -0.1)
legend("topright", title="Legend", fill = c("deeppink2", "deepskyblue3", "green3", "darkorange1"),
       legend = as.character(c("Tumor", "Normal", "Male", "Female")))

#*ArrayQualityMetrics
arrayQualityMetrics(expressionset=ExpressionSet(assayData=dat_norm_E, phenoData=annotated_targets),
                    outdir=file.path("./", "QC_norm_E"), 
                    force=TRUE, intgroup="Group")
dev.off() #exit folder

## Filtering ----
#*check normalized data
head(dat_norm)
#*remove control probes ($genes$ControlType==1)
control <- dat_norm$genes$ControlType==1
#*remove circRNAs with no name
NoName <- is.na(dat_norm$genes$circRNA)
#*remove Negative040 samples
empty <- dat_norm$genes$circRNA==""
#*check discarded probes (could be duplicated -> overestimation)
sum(control, na.rm=TRUE)
sum(NoName, na.rm=TRUE)
sum(empty, na.rm=TRUE)
#*remaining probes
dat_filter <- dat_norm[!control & !NoName & !empty, ]
#*check dataset sizes
dim(dat_norm)
dim(dat_filter)

## Filtered data quality control ----
#*expression data
dat_filter_E <- dat_filter$E

#*plotDensities 
plotDensities(dat_filter, legend="topright", main="Filtered data",
              Group=targets$Group)

#*boxplot - check signal intensities
#//las - twisted lables, range - whisker length 
#//0-all data
boxplot(dat_filter_E, range=0, las=3, ylab="log2 intensity", main="Filtered data", 
        col=my_colors)

#*PCA
plotCustomPCA(dat_filter_E, labels=targets$Name, groups=targets$Group,
              title="Filtered data", colors=c("green", "yellow", "red", "blue"),
              size=4)

#*ArrayQualityMetrics
arrayQualityMetrics(expressionset=ExpressionSet(assayData=dat_filter_E, phenoData=annotated_targets),
                    outdir=file.path("./", "QC_final_E"), 
                    force=TRUE, intgroup="Group")
dev.off() #exit folder

#*dendrogram
plotCustomDendrogram(input_data = dat_filter_E, title = "Filtered data", names = targets$Name,
                     groups = targets$Group, colors = c(Tumor="deeppink2", Normal="deepskyblue3"), 
                     groups_2 = targets$Sex, colors_2 = c(M="green3", F="darkorange1"),
                     bar_label = "Sex", hanging = -0.08)
legend("topright", title="Legend", fill = c("deeppink2", "deepskyblue3", "green3", "darkorange1"),
       legend = as.character(c("Tumor", "Normal", "Male", "Female")))

## Linear model ----
#*design matrix
design_Group <- model.matrix(~Pair+Group) #come from factors for phenoData
design_Group

#*lm fit (limma guide 8.1)
fit_Group <- lmFit(dat_filter, design_Group)

#*moderate standard errors of the estimated log-fold changes
fit2_Group <- eBayes(fit_Group)


## Find differentially expressed genes ----
#//classify results as up, down or non-significant
#*group - normal vs tumor
results_Group <- decideTests(fit2_Group)
summary(results_Group)
topTable(fit2_Group, coef = "GroupNormal") #coef - which comparison are we interested in

#make new variable
#//n=Inf for all probes that fit given parameters
all_fit2_Group <- topTable(fit2_Group, n=Inf, coef = "GroupNormal")
#// n= number of significant probes (from summary(results_Group))
top_fit2_Group_sig <- topTable(fit2_Group, n=55, coef = "GroupNormal")
#//lfc =minimum absolute log2-fold-change required
top_fit2_Group <- topTable(fit2_Group, lfc = log(1.5,2), n=Inf, coef = "GroupNormal")
#//lfc + p.value
top_fit2_Group_FC_sig <- topTable(fit2_Group, lfc = log(1.5,2), p.value = 0.05,
                                  n=Inf, coef = "GroupNormal")

#annotation for selected probes
all_fit2_Group_annotation <- match_df(probeAnnotation, all_fit2_Group, on="circRNA")
top_fit2_Group_sig_annotation <- match_df(probeAnnotation, top_fit2_Group_sig, on="circRNA")
top_fit2_Group_annotation <- match_df(probeAnnotation, top_fit2_Group, on="circRNA")
top_fit2_Group_FC_sig_annotation <- match_df(probeAnnotation, top_fit2_Group_FC_sig, on="circRNA")

#check filtering (first number must be the same)
dim(all_fit2_Group)
dim(all_fit2_Group_annotation)

dim(top_fit2_Group_sig)
dim(top_fit2_Group_sig_annotation)

dim(top_fit2_Group)
dim(top_fit2_Group_annotation)

dim(top_fit2_Group_FC_sig)
dim(top_fit2_Group_FC_sig_annotation)

#order by circRNA name
all_fit2_Group_ordered <- all_fit2_Group[order(all_fit2_Group$circRNA),]
all_fit2_Group_annotation_ordered <- all_fit2_Group_annotation[order(all_fit2_Group_annotation$circRNA),]

top_fit2_Group_sig_ordered <- top_fit2_Group_sig[order(top_fit2_Group_sig$circRNA),]
top_fit2_Group_sig_annotation_ordered <- top_fit2_Group_sig_annotation[order(top_fit2_Group_sig_annotation$circRNA),]

top_fit2_Group_ordered <- top_fit2_Group[order(top_fit2_Group$circRNA),]
top_fit2_Group_annotation_ordered <- top_fit2_Group_annotation[order(top_fit2_Group_annotation$circRNA),]

top_fit2_Group_FC_sig_ordered <- top_fit2_Group_FC_sig[order(top_fit2_Group_FC_sig$circRNA),]
top_fit2_Group_FC_sig_annotation_ordered <- top_fit2_Group_FC_sig_annotation[order(top_fit2_Group_FC_sig_annotation$circRNA),]

#select data for merging and exporting
all_final_results_Group <- cbind(all_fit2_Group_ordered[,4:4], all_fit2_Group_annotation_ordered[,2:3],
                                 all_fit2_Group_ordered[,7:11], all_fit2_Group_annotation_ordered[,4:11])

top_final_results_Group_sig <- cbind(top_fit2_Group_sig_ordered[,4:4], top_fit2_Group_sig_annotation_ordered[,2:3],
                                 top_fit2_Group_sig_ordered[,7:11], top_fit2_Group_sig_annotation_ordered[,4:11])

top_final_results_Group <- cbind(top_fit2_Group_ordered[,4:4], top_fit2_Group_annotation_ordered[,2:3],
                                 top_fit2_Group_ordered[,7:11], top_fit2_Group_annotation_ordered[,4:11])

top_final_results_Group_FC_sig <- cbind(top_fit2_Group_FC_sig_ordered[,4:4], top_fit2_Group_FC_sig_annotation_ordered[,2:3],
                                 top_fit2_Group_FC_sig_ordered[,7:11], top_fit2_Group_FC_sig_annotation_ordered[,4:11])

#rename first column
colnames(all_final_results_Group) <- c("ProbeName", colnames(all_final_results_Group)[2:16])
colnames(top_final_results_Group_sig) <- c("ProbeName", colnames(top_final_results_Group_sig)[2:16])
colnames(top_final_results_Group) <- c("ProbeName", colnames(top_final_results_Group)[2:16])
colnames(top_final_results_Group_FC_sig) <- c("ProbeName", colnames(top_final_results_Group_FC_sig)[2:16])

#export to excel
dir.create("Results", recursive = TRUE, showWarnings = FALSE)
write.xlsx(all_final_results_Group, file = "./Results/Tumor_vs_Normal_all.xlsx", row.names = FALSE)
write.xlsx(top_final_results_Group_sig, file = "./Results/Tumor_vs_Normal_significant.xlsx", row.names = FALSE)
write.xlsx(top_final_results_Group, file = "./Results/Tumor_vs_Normal_foldChange.xlsx", row.names = FALSE)
write.xlsx(top_final_results_Group_FC_sig, file = "./Results/Tumor_vs_Normal_foldChange_significant.xlsx", row.names = FALSE)

## Visualizations ----
#volcano plot
EnhancedVolcano(all_final_results_Group, #topTable data
                lab = all_final_results_Group$circRNA, #labels for data points
                x="logFC", #column with log2FC
                y="adj.P.Val", #column with p-values "P.Value" or "adj.P.Val"
                title="Differentially expressed circRNAs",
                subtitle="GSE101684 NSCLC",
                xlim=c(-2, 2),
                ylim=c(0, 2.5),
                # legendLabSize = 0,
                # legendIconSize = 0,
                pCutoff = 0.05, #horizontal line
                FCcutoff = log(1.5, 2), #vertical lines
                cutoffLineType="solid", #line type
                shape = 15, #label shape
                caption=c(), #default is paste0("Total = ", nrow(toptable), " variables")
                #label colors (default c("grey30", "forestgreen", "royalblue", "red2"))
                #< abs(FCcutoff) && > pCutoff, > abs(FCcutoff), < pCutoff, > abs(FCcutoff) && < pCutoff
                col=c("grey30", "grey30", "grey30", "red2"),
                border = "full", #full frame around graph ("full) or just x and y axis ("partial")
                gridlines.minor=FALSE #remove intermediate lines
                )

#heatmap - if error occurs, first try run dev.off() in console
sampleDend <- samplesDendrogram(input_data = dat_filter$E)
heatmap.2(dat_filter$E,
          trace="none",
          Rowv = TRUE, #reordering for dendrogram (logical or dendrogram)
          Colv = sampleDend, #reordering for dendrogram (logical or dendrogram)
          scale="row", #color scaling by each probe
          symbreaks = FALSE,
          symkey = FALSE,
          dendrogram = "both", #for both dendrograms (by sample and by probe)
          cexCol = 2, #sample label size
          main="circRNA", #main title
          col=greenred(75), #colors
          keysize = 0.75, #legend size
          #heatmap is 2x2 matrix (colour key, column dendrogram, row dendrogram, image plot)
          #ColSideColors - 2x3 matrix; RowSideColors - 3x2 matrix
          #lmat, lhei, lwid: position matrix, column height, column width
          lhei = c(1, 4, 0.005),
          #border sizes (bottom - column labels, right - row labels)
          margins = c(8,6), 
          ColSideColors = c(rep("red", 4), rep("cyan", 4)), #sample line colors
          density.info = "histogram") #legend (histogram, density or nothing)
