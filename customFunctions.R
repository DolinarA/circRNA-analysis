# Custom functions for circRNA analysis

# plotCustomPCA

# input_data=data$E ali log2(data$E)
# groups=targets$PCA_group
# labels=targets$Name
# title="Raw/Normalized/Filtered data"

plotCustomPCA <- function(input_data, labels, groups, title, 
                          colors=c("green", "yellow", "red", "blue"), 
                          scale=FALSE, size=1.5){
  data <- prcomp(t(input_data), scale=scale) # expression matrix
  
  # plot adjustments
  dataDf <- data.frame(data$x)
  Group <- groups # for legend
  percentVar <- round(data$sdev^2/sum(data$sdev^2)*100, 1) # percent of variance
  
  # main plot
  p1 <- ggplot(dataDf, aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color="gray70") +
    geom_vline(xintercept = 0, color="gray70") +
    geom_point(aes(color=Group), alpha=0.65, size=3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5, max(data$x[,1])+5)) +
    scale_fill_discrete(name="Group")
  
  # avoid label superposition
  p1+geom_text_repel(aes(y=PC2+0.25, label=labels), segment.size = 0.25, size=size) + 
    labs(x=c(paste("PC1", percentVar[1], "%")), y=c(paste("PC2", percentVar[2], "%"))) +
    ggtitle(paste("Principal Component Analysis for: ", title, sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colors)
}


# plotCustomDendrogram

# input_data=data$E or similar
# method="spearman" (default), "pearson" or "kendall"

plotCustomDendrogram <- function(input_data, method="spearman", title, names,
                                 groups, colors, groups_2, colors_2, 
                                 bar_label, hanging=0.1, label_size=1){
  
  # calculate the distances
  d <- as.dist(1-cor(input_data, method=method))
  
  # sample clustering (hclust method)
  sampleClustering <- hclust(d)
  
  # dendrogram
  sampleDendrogram <- as.dendrogram(sampleClustering, hang=hanging)
  
  # sample labels
  labels(sampleDendrogram) <- names[order.dendrogram(sampleDendrogram)]
  
  #add colors to dendrogram (for sample names) c(S="deeppink2", C="deepskyblue3")
  #groups = targets$Group
  labels_colors(sampleDendrogram) <- colors[groups][order.dendrogram(sampleDendrogram)]
  #set label size
  dend <- set(sampleDendrogram, "labels_cex", label_size) #label size
  #draw graph
  plot(dend, main=title) #"All samples" or similar
  #add colored bar (second level of differentiating samples) c(M="green3", F="darkorange1")
  #groups_2 = targets$Sex
  my_colors <- colors_2[groups_2][order.dendrogram(sampleDendrogram)]
  colored_bars(colors = my_colors, dend = dend, rowLabels = bar_label, sort_by_labels_order=FALSE)
}

#sampleDendrogram (for heatmap)
samplesDendrogram <- function(input_data, method="spearman", hanging=-1){
  d <- as.dist(1-cor(input_data, method=method)) #input_data: data$E or similar
  #sample clustering using hclust
  sampleClustering <- hclust(d)
  #dendrogram
  sampleDendrogram <- as.dendrogram(sampleClustering, hang=hanging)
}

#plotCustomRLE
plotCustomRLE <- function(input_data, colors, title,
                          x_label="Sample", y_label="log2 expression deviation"){
  dat_log2_E <- input_data #data$E or log2(data$E)
  #calculate the median for each probe (row)
  dat_log2_E_rowMedians <- rowMedians(as.matrix(dat_log2_E))
  #make discrepancy matrix (calculate the difference between expression level and median)
  RLE_raw_log2 <- sweep(dat_log2_E, 1, dat_log2_E_rowMedians)
  #convert to data.frame
  RLE_raw_log2_df <- as.data.frame(RLE_raw_log2)
  #shrink the matix
  RLE_raw_log2_df_gathered <- gather(RLE_raw_log2_df, sample_codes, log2_expression_deviation)
  #graph
  #warning: removed rows containing non-finite values - result of NA values being excluded in the plot
  ggplot(RLE_raw_log2_df_gathered, aes(sample_codes, log2_expression_deviation)) + 
    geom_boxplot(outlier.shape = NA, fill=colors) +  #fill=c(rep("yellow", 8), rep("chartreuse3", 12))
    ylim(c(-2, 2)) + 
    theme(axis.text.x = element_text(colour = "aquamarine4", angle = 90, size = 8,
                                     hjust = 0.1, face = "bold"),
          plot.title = element_text(hjust = 0.5)) +
    labs(title=paste("Relative log2 expression -", title, sep=" "), x=x_label, y=y_label) #title="raw/background corrected data"
}

#plotCustomPVCA
plotCustomPVCA <- function(input_data, pheno_data, batch_factors, threshold = 0.6, x_lab="Effects",
                           y_lab="Weighted average proportion variance",
                           graph_title="Batch effect estimation\n(Principal variance component analysis method)"){
  #convert data to eSet (input_data=data$E, pheno_data=annotated_targets)
  eset_dat_norm <- ExpressionSet(assayData = input_data, phenoData = pheno_data)
  pct_threshold <- threshold #select threshold
  batch.factors <- batch_factors #c("Group", "Sex", "C9", "Age")
  pvcaObj <- pvcaBatchAssess(eset_dat_norm, batch.factors, pct_threshold)
  
  #set plot margins
  #bottom, left, top, and right; the default is c(5.1, 4.1, 4.1, 2.1).
  par(mar=c(7, 4.1, 4.1, 2.1)) 
  bp <- barplot(pvcaObj$dat, ylab=y_lab,
                ylim=c(0, 1.0),
                col="mediumorchid", las=2, 
                main=graph_title)
  #x-axis labels
  axis(1, at=bp, labels = pvcaObj$label, cex.axis=0.9, las=2)
  title(xlab=x_lab, line=5.5) #line changes the distance between the title and the axis
  #add % of variance
  dat_values=pvcaObj$dat
  new_values=round(dat_values, 3)
  text(bp, pvcaObj$dat, labels = new_values, pos=3, cex=0.9)
}