options(stringsAsFactors = F)
library(grid)
library(gridExtra)
library(plyr)
library(ggplot2)
library(sva)


# Analisi espressione campioni ===============


f  <- list.files(path="expression",
                 recursive=F,
                 pattern=".gene.counts"
                 ,full.names=T)


a=read.delim2(f[1], header=T, skip=1)


map=a[,c("Geneid", "Length","gene_name", "gene_type")]


map2=map

for(i in f){
  
  b=read.delim2(i, header=T, skip=1)
  print(dim(b))
  map2$new=b[,ncol(b)]
  map2$new=as.numeric(map2$new)
  colnames(map2)[ncol(map2)]=gsub("\\.", "-", gsub(".Aligned.sortedByCoord.out.bam", "", colnames(b)[ncol(b)]))
}


rownames(map2)=map2$Geneid

saveRDS(map2, "Rdata/expression.rds")



# Batch =========
c=readRDS("Rdata/expression.rds")
c=subset(c, gene_type=="protein_coding")

info=read.csv("Tables/Scheme.csv")

batch=info$rep[match(colnames(c)[5:12], info$barcode)]

c2=ComBat_seq(as.matrix(c[,5:12]), batch=batch)
c2=cbind.data.frame(c[, 1:4], c2)


saveRDS(c2, "Rdata/Gene_expression_batch_protein_coding.rds")


# Init ####
library(tidyverse)

get_DESEQ = function(counts, samples, cond, ctrl, tlen, cores=1){
  
  # Perform DESeq analysis
  require(DESeq2)
  require(BiocParallel)
  
  tmp= floor(as.matrix(counts[,samples]))
  colData <- data.frame(condition=factor(cond))
  colData$condition <- relevel(colData$condition, ref = ctrl) 
  dds <- DESeqDataSetFromMatrix(tmp, colData, formula(~condition))
  mcols(dds)$basepairs=tlen
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds, parallel = T, BPPARAM = MulticoreParam(workers = cores))
  dds
}

get_DESEQ_norm = function(counts, samples, cond, ctrl, tlen, cores=1){
  
  # Perform DESeq analysis
  require(DESeq2)
  require(BiocParallel)
  
  tmp= floor(as.matrix(counts[,samples]))
  colData <- data.frame(condition=factor(cond))
  colData$condition <- relevel(colData$condition, ref = ctrl) 
  dds <- DESeqDataSetFromMatrix(tmp, colData, formula(~condition))
  dds
}

plotRNAVolcanos <- function(de, lfcTh=1, pvTh=0.05, top=14)
{
  require(ggplot2)
  require(RColorBrewer)
  require(ggrepel)
  
  pal <- brewer.pal(9, "Set1")
  
  fcIdx <- grep("logFC|log2FoldChange", colnames(de), value = T)
  pvIdx <- grep("adj.P.Val|FDR|padj", colnames(de), value = T)
  tmp   <- de[,c(fcIdx,pvIdx)]
  colnames(tmp) <- c("lfc","padj")
  tmp$status <- "none"
  tmp$status[which(tmp$lfc>=lfcTh & tmp$padj<=pvTh)] <- "Up-regulated"
  tmp$status[which(tmp$lfc<=(-lfcTh) & tmp$padj<=pvTh)] <- "Down-regulated"
  tmp$status <- factor(tmp$status)
  
  topUP <- tmp[tmp[,'status']=="Up-regulated",]
  topUP <- rownames(topUP[1:top,])
  topDW <- tmp[tmp[,'status']=="Down-regulated",]
  topDW <- rownames(topDW[1:top,])
  
  tmp$lab <- rownames(tmp)
  tmp$lab[-which(rownames(tmp)%in%c(topDW, topUP))] <- ""
  
  p <- ggplot(tmp, aes(x=lfc, y=-log10(padj),col=status, label=lab)) + geom_point(size=2) + theme_bw() +
    geom_hline(yintercept = -log10(pvTh), linetype = 'dashed') +
    geom_vline(xintercept = lfcTh, linetype = 'dashed') +
    geom_vline(xintercept = -lfcTh, linetype = 'dashed') +
    xlab("logFC") + ylab("adjusted P-value") +
    theme_bw() + theme(panel.grid = element_blank(),legend.position = "right",
                       aspect.ratio = 0.8, plot.title = element_text(face="bold", hjust = 0.5)) +
    scale_color_manual(values =c("grey","Up-regulated"="red" )) +
    geom_text_repel()
  
  return(p)
}


# DE======

c=readRDS("Rdata/Gene_expression_batch_protein_coding.rds")

info=read.csv("Tables/Scheme.csv")


rownames(c) = paste0(c$gene_name,'_', rownames(c))

samples=colnames(c)[grep("Rotolo", colnames(c))]

group=info$group[match(samples, info$barcode)]

deg = get_DESEQ(c
                , samples=samples
                , cond = group
                , ctrl='group1'
                , c$Length)

res <- results(deg, name="condition_group2_vs_group1")
res <- res[order(res$pvalue),]


res$gene_name=c$gene_name[match(rownames(res), rownames(c))]
res$gene_id=c$Geneid[match(rownames(res), rownames(c))]

write.csv(res, "Tables/Differential_gene_expression_batch_protein_coding.csv", row.names = F)

# Figures =============


res2=as.data.frame(res)
res2=subset(res2, !is.na(padj))

rownames(res2)=res2$gene_name

pdf(file='Figures/Volcano_batch_protein_coding.pdf', height = unit(10,'cm'), width = unit(11,'cm'), useDingbats = F )
plotRNAVolcanos(res2)
dev.off()


sig=subset(as.data.frame(res),padj<=0.05 & abs(log2FoldChange)>1 )

nor=as.data.frame(assay(vsd))
nor$gene_name=c$gene_name[match(rownames(nor), rownames(c))]
nor$gene_id=c$Geneid[match(rownames(nor), rownames(c))]

nor=subset(nor, gene_name%in%sig$gene_name)

rownames(nor)=nor$gene_name
nor$gene_id<-NULL
nor$gene_name<-NULL

nor=as.matrix(nor)


names(group)=samples
m=nor
annotCol =group
fig_h = unit(4,'cm')
fig_w = unit(4,'cm')
retHm    = F

  
require(ComplexHeatmap)
require(circlize)
require(RColorBrewer)
  
base_mean <- rowMeans(m)
m_scaled <- t(apply(m, 1, scale))
colnames(m_scaled) <- colnames(m)
  
bPalette <- colorRampPalette(brewer.pal(9, "Reds"))
myPalette <- c("blue","black","red")
ramp <- colorRamp2(c(-2, 0, 2), myPalette)
  
rep=info$rep
names(rep)=samples

ha_column <- HeatmapAnnotation(group= annotCol, rep=rep, 
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8))
                                     ,annotation_height = unit(4, "mm"),
                               col = list(group=c("group1"="purple", "group2"="forestgreen"), rep=c("rep1"="magenta", "rep2"="cyan", "rep3"="yellow", "rep4"="orange")),
                               border = T)

hm <- Heatmap(m_scaled,
                col = ramp,
                show_row_dend = T,
                row_names_side = "left",
                row_names_gp = gpar(fontsize=8),
                column_names_gp = gpar(fontsize=8),
                column_title_gp = gpar(fontsize=10, fontface="bold"),
                heatmap_legend_param = list(title = "row Z-score",
                                            title_gp = gpar(fontsize=8),
                                            title_position = "topcenter",
                                            # legend_width  = unit(4, "cm"),
                                            # legend_height = unit(0.5, "mm"),
                                            values_gp     = gpar(fontsize=8),
                                            legend_direction = "vertical")
                , top_annotation = ha_column
                )
  
bmscale <- summary(base_mean)
bmramp <- colorRamp2(c(bmscale[1],bmscale[3],bmscale[5]), bPalette(3))
bmh <- Heatmap(base_mean
                 , name = "Mean Expression"
                 , column_names_gp = gpar(fontsize=8)
                 , show_row_names = FALSE
                 , width = unit(3, "mm")
                 , col = bmramp
                 , heatmap_legend_param = list(title = "Base Mean",title_gp = gpar(fontsize=8)))
  
hmOut <- hm + bmh


pdf(file='Figures/Heatmap_padj_0_05_FC_1_batch_protein_coding.pdf', height = unit(4.5,'cm'), width = unit(5,'cm'), useDingbats = F )
hmOut
dev.off()
