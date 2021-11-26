
##no renv, chunks for knitr


## ---- packages


options(connectionObserver = NULL) # https://github.com/rstudio/rstudio/issues/9219

library(mouse4302.db)
library(pd.mouse430.2)
library(org.Mm.eg.db) # issues with loading

library(oligo)
library(arrayQualityMetrics)

library(limma)
library(edgeR)
library(multiMiR)

library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(knitr)
library(xtable)
library(rlist)

library(ComplexHeatmap)
library(circlize)

library(genefilter)
#library(geneplotter)
library(topGO)
library(Rgraphviz)

library(reactome.db)
library(ReactomePA)
library(DOSE)
library(enrichplot)

#format
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
fontsize2 <- theme(axis.text=element_text(size=20), axis.title=element_text(size=20))
theme_pca<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    strip.background=element_blank(),axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
    legend.text=element_text(size=10) )


## ---- prep_environment
projdir="/Users/agata.smialowska/NBISproj/1186/analysis/report"


scriptdir=file.path(projdir,"scripts")
func_pth=file.path(scriptdir,"sel_mir_fnctions.R")
source(func_pth)

resdir=file.path(projdir,"results_prel_v.1.0_15vi2021")
dir.create(resdir)



## ---- mRNA_data

data_cel="/Users/agata.smialowska/NBISproj/1186/GSE12049_RAW"

res_mRNA_dir=file.path(resdir,"mRNA_GSE12049")
dir.create(res_mRNA_dir)

Array.Data.File=c("GSM304403.CEL","GSM304404.CEL","GSM304405.CEL","GSM304406.CEL","GSM304407.CEL","GSM304408.CEL")
Factor.Value.phenotype=c(rep("WT",3),rep("Mut",3))

samples=cbind(Array.Data.File,Factor.Value.phenotype)
samples=as.data.frame(samples)
rownames(samples)=samples$Array.Data.File


# > samples
#               Array.Data.File Factor.Value.phenotype
# GSM304403.CEL   GSM304403.CEL                     wt
# GSM304404.CEL   GSM304404.CEL                     wt
# GSM304405.CEL   GSM304405.CEL                     wt
# GSM304406.CEL   GSM304406.CEL                     kn
# GSM304407.CEL   GSM304407.CEL                     kn
# GSM304408.CEL   GSM304408.CEL                     kn


#read in data
raw_data = oligo::read.celfiles(filenames = file.path(data_cel, samples$Array.Data.File))

#add annotation
pheno = AnnotatedDataFrame(data= samples)
phenoData(raw_data) = pheno


# > raw_data
# ExpressionFeatureSet (storageMode: lockedEnvironment)
# assayData: 1004004 features, 6 samples 
#   element names: exprs 
# protocolData
#   rowNames: GSM304403.CEL GSM304404.CEL ... GSM304408.CEL (6 total)
#   varLabels: exprs dates
#   varMetadata: labelDescription channel
# phenoData
#   rowNames: GSM304403.CEL GSM304404.CEL ... GSM304408.CEL (6 total)
#   varLabels: Array.Data.File Factor.Value.phenotype
#   varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: pd.mouse430.2 


## ---- mRNA_raw_exploration_QC

exp_raw = log2(Biobase::exprs(raw_data))
PCA_raw = prcomp(t(exp_raw), scale. = FALSE)

percentVar = round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio = sqrt(percentVar[2] / percentVar[1])

data_pca = data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                    Strain = pData(raw_data)$Factor.Value.phenotype)

pca_raw_mRNA=ggplot(data_pca, aes(PC1, PC2)) +
      geom_point(aes(colour = Strain)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4")) + fontsize + theme_pca


intensities_raw=oligo::boxplot(exprs(raw_data), which="pm", target="core",
               main = "Boxplot of log2-intensitites for the raw data")


qcdir=file.path(res_mRNA_dir,"arrayQC")
dir.create(qcdir)

arrayQualityMetrics(expressionset = raw_data,
    outdir = qcdir,
    force = TRUE, do.logtransform = TRUE,
    intgroup = "Factor.Value.phenotype")

## ---- mRNA_normalisation

#RLE before normalisation
data_rma = oligo::rma(raw_data, normalize = FALSE)

row_medians_assayData = Biobase::rowMedians(as.matrix(Biobase::exprs(data_rma)))

RLE_data = sweep(Biobase::exprs(data_rma), 1, row_medians_assayData)
RLE_data = as.data.frame(RLE_data)
RLE_data_gathered = tidyr::gather(RLE_data, sample, log2_expression_deviation)

box_int_raw=ggplot2::ggplot(RLE_data_gathered, aes(sample, log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-1, 1.5)) + 
  theme(axis.text.x = element_text(angle = 60, size = 8, hjust = 1 ,face = "bold"))


## normalisation
data_rma_norm = oligo::rma(raw_data)

row_medians_assayData = Biobase::rowMedians(as.matrix(Biobase::exprs(data_rma_norm)))

RLE_data_norm = sweep(Biobase::exprs(data_rma_norm), 1, row_medians_assayData)
RLE_data_norm = as.data.frame(RLE_data_norm)
RLE_data_gathered_norm = tidyr::gather(RLE_data_norm, sample, log2_expression_deviation)

box_int_rma=ggplot2::ggplot(RLE_data_gathered_norm, aes(sample,log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-1, 1.5)) +
  theme(axis.text.x = element_text(angle = 60, size = 8, hjust = 1 ,face = "bold"))

intensities_rma=oligo::boxplot(exprs(data_rma_norm), which="pm", target="core",
               main = "Boxplot of log2-intensitites for the rma-normalised data")


## PCA
exp_norm = Biobase::exprs(data_rma_norm)
PCA_norm = prcomp(t(exp_norm), scale. = FALSE)

percentVar = round(100*PCA_norm$sdev^2/sum(PCA_norm$sdev^2),1)
sd_ratio = sqrt(percentVar[2] / percentVar[1])

data_pca = data.frame(PC1 = PCA_norm$x[,1], PC2 = PCA_norm$x[,2],
                   Strain = pData(raw_data)$Factor.Value.phenotype)

pca_norm_mRNA=ggplot(data_pca, aes(PC1, PC2)) +
      geom_point(aes(colour = Strain), size=4) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4")) + fontsize + theme_pca


## ---- box_int_raw
box_int_raw + theme_bw()

## ---- box_int_rma
box_int_rma + theme_bw()

## ---- pca_microarray
pca_norm_mRNA


## ---- mRNA_probe_filtering

# filtering low intensity probes
medians = rowMedians(Biobase::exprs(data_rma_norm))

man_threshold = 4

hist_filt = hist(medians, 100, col = "cornsilk", freq = FALSE, 
            main = "Histogram of the median intensities",
            border = "antiquewhite4",
            xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)

samples_cutoff=3

# find indices which should be kept
idx_man_threshold = apply(Biobase::exprs(data_rma_norm), 1,
                           function(x){
                          sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

# idx_man_threshold
# FALSE  TRUE 
# 18106 26995 

probes_low_intensity=length(idx_man_threshold[idx_man_threshold==FALSE])

data_manfilter_keep = subset(data_rma_norm, idx_man_threshold)


# filtering probes which associate to more than one gene

# the following did not work recently due to issues with loading od org.Mm.eg.db (solved though, see above)
annotation = AnnotationDbi::select(mouse4302.db,
                                  keys = (featureNames(data_manfilter_keep)),
                                  columns = c("SYMBOL", "GENENAME","ENSEMBL"),
                                  keytype = "PROBEID")


# filter annotations
annotation = subset(annotation, !is.na(SYMBOL))
anno_grouped = group_by(annotation, PROBEID)
anno_summarized = dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

# > nrow(anno_summarized)
# [1] 25247

# remove probes with more than one gene annotation
anno_filtered_out = filter(anno_summarized, no_of_matches > 1) 

# > nrow(anno_filtered_out)
# [1] 882

probes_multi_annot=nrow(anno_filtered_out)

ids_to_exlude = (featureNames(data_manfilter_keep) %in% anno_filtered_out$PROBEID)

table(ids_to_exlude)
# ids_to_exlude
# FALSE  TRUE 
# 26113   882 


data_final = subset(data_manfilter_keep, !ids_to_exlude)

#> nrow(data_final)
# Features 
#    26113 

probes_retained=nrow(data_final)


fData(data_final)$PROBEID = rownames(fData(data_final))
annot_final=merge(fData(data_final),annotation,by="PROBEID",all.x=TRUE)
annot_final=annot_final[!duplicated(annot_final$PROBEID),]

#subset data to match filtered annot
annot_final_data=left_join(fData(data_final), annot_final)
rownames(annot_final_data)=annot_final_data$PROBEID

fData(data_final)=annot_final_data

# > validObject(data_final)
# [1] TRUE




## ---- mRNA_DE_show

# using limma

group=as.factor(samples$Factor.Value.phenotype)

mod.mat=model.matrix(~ 0 + group)
colnames(mod.mat)=c("Mut","WT")
rownames(mod.mat)=rownames(samples)



mod.mat

# > mod.mat
#               kn wt
# GSM304403.CEL  0  1
# GSM304404.CEL  0  1
# GSM304405.CEL  0  1
# GSM304406.CEL  1  0
# GSM304407.CEL  1  0
# GSM304408.CEL  1  0
# attr(,"assign")
# [1] 1 1
# attr(,"contrasts")
# attr(,"contrasts")$group
# [1] "contr.treatment"

contrast.mat = makeContrasts(Mut-WT, levels = mod.mat)

contrast.mat

#       Contrasts
# Levels kn - wt
#     kn       1
#     wt      -1

## ---- mRNA_DE


fit = eBayes(contrasts.fit(lmFit(data_final,
                                design = mod.mat),
                                contrast.mat))


results_mRNA = topTable(fit, number = nrow(fit))

fname="mRNA-GSE12049-MutVsWT.limma.tsv"
fpth=file.path(res_mRNA_dir,fname)
write.table(results_mRNA,fpth,sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


hist_pvals_marray=hist(results_mRNA$P.Value, col = "cornsilk",
     main = "p values for comparison KD vs WT", xlab = "p-values")



#results_mRNA_sig=results_mRNA[results_mRNA$adj.P.Val<0.05,]
#> nrow(results_mRNA_sig)
# [1] 238


#results_mRNA_sig=results_mRNA[results_mRNA$adj.P.Val<0.1,]
# > nrow(results_mRNA_sig)
# [1] 840

############
# selection for miRNA targets overlap

# significance
# arbitrary FDR
#FDR_CO=0.05
FDR_CO=0.1

#FDR distribution

# > summary(results_mRNA$adj.P.Val)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0002013 0.3605238 0.6217462 0.5921558 0.8403753 0.9997648 


#FDR_CO=quantile(results_mRNA$adj.P.Val, probs = c(0,0.01,0.05,0.1,0.25,0.5), na.rm = FALSE,names = TRUE, type=1)[3]

#           0%           1%           5%          10%          25%          50% 
# 0.0002013071 0.0503282484 0.1242627888 0.1933882516 0.3605238216 0.6217461754 



#effect size
#lfc_CO=0.58 #50% change, i.e. log2(1.5)=0.5849625

lfc_CO=0 #any change

#lfc_CO=0.58 #change 1.5x


lfc_CO_lin=2^lfc_CO


df_volcano=as.data.frame(cbind(results_mRNA$logFC, results_mRNA$adj.P.Val))
colnames(df_volcano)=c("log2FoldChange","FDR")

p_volcano = ggplot(data=df_volcano, aes(x=log2FoldChange, y=-log10(FDR))) + geom_point() + theme_minimal() +
  geom_vline(xintercept=c(-lfc_CO, lfc_CO), col="red") +
  geom_hline(yintercept=-log10(FDR_CO), col="red")



results_mRNA_sig=results_mRNA[results_mRNA$adj.P.Val<FDR_CO,]

mRNA_up=results_mRNA_sig[results_mRNA_sig$logFC>lfc_CO,]
mRNA_dn=results_mRNA_sig[results_mRNA_sig$logFC<(-lfc_CO),]

# FDR 0.05 lFC 0
# > nrow(mRNA_up)
# [1] 127
# > nrow(mRNA_dn)
# [1] 111


## ---- hist_FDR
hist(results_mRNA$adj.P.Val, main="", xlab="BH-adjusted p values")
abline(v=FDR_CO, col="red",lwd=3)


## ----volcano_mRNA
p_volcano



## ---- table_mRNA_res
#table for the report
results_mRNA_sig_tab=results_mRNA_sig[order(abs(results_mRNA_sig$logFC), decreasing=TRUE),][1:10,]
results_mRNA_sig_tab$logFC=as.numeric(results_mRNA_sig_tab$logFC)

caption.tab="Top 10 differentially expressed genes in Mut vs. WT, ranked by the absolute log fold change."

cols_to_print=as.vector(c("SYMBOL","GENENAME","logFC","adj.P.Val"))
rws=seq(1, (nrow(results_mRNA_sig_tab)-1), by = 2)
col=rep("\\rowcolor[gray]{0.95}", length(rws))
print(xtable(as.data.frame(results_mRNA_sig_tab[,cols_to_print]) , align=c("p{0.1cm}","c","c","c","c"), 
display=c("d","s","s","f","e"), digits=c(0,0,0,3,3),
caption=caption.tab),
include.rownames=FALSE,floating=TRUE,
table.placement="H", caption.placement = "bottom",
add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")



## ---- hm_DE_mRNA_prep
results_mRNA_sig_hm=rownames(results_mRNA_sig[order(abs(results_mRNA_sig$logFC), decreasing=TRUE),][1:30,]) #probeIDs
exp_norm_hm=exp_norm[rownames(exp_norm)%in%results_mRNA_sig_hm,]
exp_norm_hm_scaled=exp_norm_hm-rowMeans(exp_norm_hm)

#change probeIDS to gene IDS
probeIDs=rownames(exp_norm_hm_scaled)
sel=annot_final[,c(1,2)][annot_final$PROBEID%in%probeIDs,]
#this above preserves the order of elements
#make sure that the order of sel and probeids from rownames is the same
sel_reord=sel[match(probeIDs,sel$PROBEID),]
rownames(exp_norm_hm_scaled)=sel_reord$SYMBOL

colour_palette_hm=colorRamp2(c(-4, 0, 4), c("dodgerblue4", "ivory", "chocolate3")) #z-scored

annot_col=samples
annot_col=as.data.frame(annot_col[,-1])
colnames(annot_col)="Strain"

ann_colors=list(
	Strain=c(Mut="darkorange2", WT="dodgerblue4"))

hm_mRNA=ComplexHeatmap::pheatmap(exp_norm_hm_scaled ,color = colour_palette_hm,
	fontsize=6, fontsize_row=7,cellheight=7.5,
	annotation_col=annot_col,annotation_colors=ann_colors)

graphics.off()

## ---- hm_DE_mRNA
print(hm_mRNA)



## ---- plots_marray_save
graphics.off()
fname="mRNA-GSE12049-plots.pdf"
fpth=file.path(res_mRNA_dir,fname)
pdf(fpth)
pca_raw_mRNA
intensities_raw
box_int_raw
box_int_rma
intensities_rma
pca_norm_mRNA +  ggtitle("PCA plot of the log-transformed \n rma normalised mRNA expression")
hist_filt
hist_pvals_marray
hm_mRNA
dev.off()


##########
# miRNA DE

## ---- miRNA_data

datadir="/Users/agata.smialowska/NBISproj/1186/gitlab/muscle_mirna/src"

res_miRNA_dir=file.path(resdir,"miRNA_up072")
dir.create(res_miRNA_dir)


meta="metadata_up_072.txt"
counts="up_072_count.tsv"

meta_pth=file.path(datadir,meta)
counts_pth=file.path(datadir,counts)

meta.data=read.delim(meta_pth,sep="\t")
meta.data$condition=c(rep("Mut",5),rep("WT",5))

mir.counts=read.delim(counts_pth,sep="\t")
rownames(mir.counts)=mir.counts$miR
mir.counts=mir.counts[,-1]

# > dim(mir.counts)
# [1] 1907   10



## ---- miRNA_DGE_exploration

group=as.factor(meta.data$condition)

mir.dge = DGEList(counts=mir.counts, group=group)

keep = filterByExpr(mir.dge)

summary(keep)
#    Mode   FALSE    TRUE 
# logical    1540     367 

nMiR_exprs=length(keep[keep=="TRUE"])


mir.dge = mir.dge[keep, , keep.lib.sizes=FALSE]

mir.dge = calcNormFactors(mir.dge)

# $samples
#           group lib.size norm.factors
# up_072_1    mut  3567130    1.6791367
# up_072_2    mut  4810590    1.5931048
# up_072_3    mut  5925672    1.3819703
# up_072_4    mut  3141630    1.4420395
# up_072_5    mut  4675535    1.4032843
# up_072_6    ref 10168859    0.6246392
# up_072_7    ref  9475429    0.7410713
# up_072_8    ref  4433827    0.6720739
# up_072_9    ref  7822736    0.5967125
# up_072_10   ref  6716639    0.7200734


mir.dge = estimateDisp(mir.dge)

# mir.dge$common.dispersion
# [1] 0.07327371

#BCV
# > sqrt(mir.dge$common.dispersion)
# [1] 0.2706912

#https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#p 23
#Typical values for the common BCV (square-root-dispersion) for datasets arising from well-controlled experiments are 0.4 for human data,0.1 for data on genetically identical model organisms or 0.01 for technical replicates.


## PCA

exp_mir <- cpm(mir.dge, normalized.lib.sizes = TRUE,log = TRUE, prior.count = 1)


PCA_mir <- prcomp(t(exp_mir), scale = FALSE, center=TRUE)

percentVar <- round(100*PCA_mir$sdev^2/sum(PCA_mir$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

data_pca <- data.frame(PC1 = PCA_mir$x[,1], PC2 = PCA_mir$x[,2],
                    Sample = meta.data$condition)


pca_mir=ggplot(data_pca, aes(PC1, PC2)) +
      geom_point(aes(colour = Sample), size=4) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4")) + fontsize + theme_pca + coord_fixed()



## ---- pca_mir
pca_mir

## ---- miRNA_DE


mir.res = exactTest(mir.dge,pair=c("WT","Mut"))

topTags(mir.res)


# Comparison of groups:  mut-ref 
#                     logFC    logCPM       PValue          FDR
# mmu-miR-544-3p  -2.720783  6.313818 1.177936e-42 4.323026e-40
# mmu-miR-3544-5p -2.162887  6.037488 8.734773e-30 1.602831e-27
# mmu-miR-1a-3p   -2.424925 18.618462 1.007885e-25 1.232980e-23
# mmu-miR-337-3p  -2.041210 11.686827 4.762982e-21 4.370036e-19
# mmu-miR-22-3p   -2.059126 10.955544 7.513769e-21 5.515107e-19



results_miRNA=as.data.frame(topTags(mir.res,n=nrow(mir.res)))
# > head(res)
#                     logFC    logCPM       PValue          FDR
# mmu-miR-544-3p   2.720783  6.313818 1.177936e-42 4.323026e-40
# mmu-miR-3544-5p  2.162887  6.037488 8.734773e-30 1.602831e-27
# mmu-miR-1a-3p    2.424925 18.618462 1.007885e-25 1.232980e-23
# mmu-miR-337-3p   2.041210 11.686827 4.762982e-21 4.370036e-19
# mmu-miR-22-3p    2.059126 10.955544 7.513769e-21 5.515107e-19


# > nrow(foo[foo$FDR<0.01,])
# [1] 149

results_miRNA$miR=rownames(results_miRNA)

results_miRNA=results_miRNA[,c(5,1:4)]

fname="miRNA-up072-MutVsWT.edgeR.tsv"
fpth=file.path(res_miRNA_dir,fname)
write.table(results_miRNA,fpth,sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


hist_pvals_mir=hist(results_miRNA$FDR, col = "cornsilk",
     main = "FDR for comparison KD vs WT", xlab = "FDR")


## subsets
FDR_COmir=0.01
lfc_COmir=0

mir_sig=results_miRNA[results_miRNA$FDR<FDR_COmir ,]

mir_up=mir_sig[mir_sig$logFC>lfc_COmir,]
mir_dn=mir_sig[mir_sig$logFC<(-lfc_COmir),]

# > nrow(mir_up)
# [1] 53
# > nrow(mir_dn)
# [1] 96


## ---- table_miRNA_res
#table for the report
mir_sig_tab=mir_sig[order(abs(mir_sig$logFC), decreasing=TRUE),][1:10,]
mir_sig_tab$logFC=as.numeric(mir_sig_tab$logFC)

caption.tab="Top 10 differentially expressed miRNAs in Mut vs. WT, ranked by the absolute log fold change."

cols_to_print=as.vector(c("miR","logFC","FDR"))
rws=seq(1, (nrow(mir_sig_tab)-1), by = 2)
col=rep("\\rowcolor[gray]{0.95}", length(rws))
print(xtable(as.data.frame(mir_sig_tab[,cols_to_print]) , align=c("p{0.1cm}","c","c","c"), 
display=c("d","s","f","e"), digits=c(0,0,3,3),
caption=caption.tab),
include.rownames=FALSE,floating=TRUE,
table.placement="H", caption.placement = "bottom",
add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")

## ---- hm_DE_miRs_prep
mir_sig_hm=rownames(mir_sig[order(abs(mir_sig$logFC), decreasing=TRUE),][1:30,])
exp_mir_hm=exp_mir[rownames(exp_mir)%in%mir_sig_hm,]
exp_mir_hm_scaled=exp_mir_hm-rowMeans(exp_mir_hm)

colour_palette_hm2=colorRamp2(c(-2, 0, 2), c("dodgerblue4", "ivory", "chocolate3")) #z-scored

annot_col=meta.data
rownames(annot_col)=meta.data$library
annot_col=as.data.frame(annot_col[,-1])
colnames(annot_col)="Strain"

ann_colors=list(
	Strain=c(Mut="darkorange2", WT="dodgerblue4"))

hm_mirs=ComplexHeatmap::pheatmap(exp_mir_hm_scaled, scale="row",color = colour_palette_hm2,
	fontsize=6, fontsize_row=7,cellheight=7.5,
	annotation_col=annot_col,annotation_colors=ann_colors)


## ---- hm_DE_miRs
print(hm_mirs)





## ---- plots_miR_save

fname="miRNA-up072-plots.pdf"
fpth=file.path(res_miRNA_dir,fname)
pdf(fpth)
pca_mir +ggtitle("PCA plot of the log-transformed miRNA expression data")
hist_pvals_mir
hm_mirs
dev.off()

###########################
##########################
## miRNA - mRNA targets


## ---- get_targets_background_not_in_the_report

#this takes some time!

# for background sets, need to be split because of query size
# all genes on the array which are targeted by miRNA in the analysis

#  length(results_mRNA$SYMBOL)
# [1] 26113
# > length(unique(results_mRNA$SYMBOL))
# [1] 14498

# > nrow(results_miRNA)
# [1] 367

## validated

all_genes_array_targeted_miRNA_valid_1=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[1:5000],
                         table   = "validated",
                         summary = TRUE,
                         use.tibble = TRUE)


all_genes_array_targeted_miRNA_valid_2=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[5001:8000],
                         table   = "validated",
                         summary = TRUE,
                         use.tibble = TRUE)

all_genes_array_targeted_miRNA_valid_3=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[8001:length(unique(results_mRNA$SYMBOL))],
                         table   = "validated",
                         summary = TRUE,
                         use.tibble = TRUE)


# predicted

#top 10 of conserved

all_genes_array_targeted_miRNA_pred_1=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[1:5000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE)

all_genes_array_targeted_miRNA_pred_2=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[5001:8000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE)

all_genes_array_targeted_miRNA_pred_3=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[8001:length(unique(results_mRNA$SYMBOL))],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE)



# top20 of all

a20_all_genes_array_targeted_miRNA_pred_1=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[1:2000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE,
                         predicted.site  = "all")

a20_all_genes_array_targeted_miRNA_pred_2=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[2001:5000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE,
                         predicted.site  = "all")

a20_all_genes_array_targeted_miRNA_pred_3=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[5001:8000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE,
                         predicted.site  = "all")


a20_all_genes_array_targeted_miRNA_pred_4=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[8001:10000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE,
                         predicted.site  = "all")

a20_all_genes_array_targeted_miRNA_pred_5=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[10001:length(unique(results_mRNA$SYMBOL))],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE,
                         predicted.site  = "all")


# top10 of all

a10_all_genes_array_targeted_miRNA_pred_1=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[1:2000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE,
                         predicted.site  = "all")

a10_all_genes_array_targeted_miRNA_pred_2=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[2001:5000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE,
                         predicted.site  = "all")

a10_all_genes_array_targeted_miRNA_pred_3=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[5001:8000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE,
                         predicted.site  = "all")


a10_all_genes_array_targeted_miRNA_pred_4=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[8001:10000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE,
                         predicted.site  = "all")

a10_all_genes_array_targeted_miRNA_pred_5=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[10001:length(unique(results_mRNA$SYMBOL))],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 10,
                         use.tibble = TRUE,
                         predicted.site  = "all")






#top 20 of conserved



c20_all_genes_array_targeted_miRNA_pred_1=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[1:5000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE)

c20_all_genes_array_targeted_miRNA_pred_2=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[5001:8000],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE)

c20_all_genes_array_targeted_miRNA_pred_3=get_multimir(org  = "mmu",
                         mirna   = rownames(results_miRNA),
                         target  = unique(results_mRNA$SYMBOL)[8001:length(unique(results_mRNA$SYMBOL))],
                         table   = "predicted",
                         summary = TRUE,
                         predicted.cutoff.type = "p",
                         predicted.cutoff      = 20,
                         use.tibble = TRUE)





#save all
targets_dir=file.path(projdir,"mRNA_GSE12049_miRNA_072_targets_multimir")
dir.create(targets_dir)

valid_targ_pth=file.path(targets_dir,"data.validated.targets.7vi2021.RData")
save(all_genes_array_targeted_miRNA_valid_1,all_genes_array_targeted_miRNA_valid_2,all_genes_array_targeted_miRNA_valid_3, file = valid_targ_pth)

pred_targ_pth=file.path(targets_dir,"data.predicted.targets.7vi2021_top20_all.RData")
save(a20_all_genes_array_targeted_miRNA_pred_1,a20_all_genes_array_targeted_miRNA_pred_2,a20_all_genes_array_targeted_miRNA_pred_3,a20_all_genes_array_targeted_miRNA_pred_4,a20_all_genes_array_targeted_miRNA_pred_5, file = pred_targ_pth)

pred_targ_pth=file.path(targets_dir,"data.predicted.targets.7vi2021_top10_conserved.RData")
save(all_genes_array_targeted_miRNA_pred_1,all_genes_array_targeted_miRNA_pred_2,all_genes_array_targeted_miRNA_pred_3, file = pred_targ_pth)

pred_targ_pth=file.path(targets_dir,"data.predicted.targets.7vi2021_top10_all.RData")
save(a10_all_genes_array_targeted_miRNA_pred_1,a10_all_genes_array_targeted_miRNA_pred_2,a10_all_genes_array_targeted_miRNA_pred_3,a10_all_genes_array_targeted_miRNA_pred_4,a10_all_genes_array_targeted_miRNA_pred_5, file = pred_targ_pth)



pred_targ_pth=file.path(targets_dir,"data.predicted.targets.7vi2021_top20_conserved.RData")
save(c20_all_genes_array_targeted_miRNA_pred_1,c20_all_genes_array_targeted_miRNA_pred_2,c20_all_genes_array_targeted_miRNA_pred_3, file = pred_targ_pth)


## ---- get_targets_background


targets_dir=file.path(projdir,"mRNA_GSE12049_miRNA_072_targets_multimir")

object_fname="data.validated.targets.7vi2021.RData"
object_fname_pth=file.path(targets_dir,object_fname)
load(object_fname_pth)


object_fname="data.predicted.targets.7vi2021_top10_all.RData"
object_fname_pth=file.path(targets_dir,object_fname)
load(object_fname_pth)

object_fname="data.predicted.targets.7vi2021_top10_conserved.RData"
object_fname_pth=file.path(targets_dir,object_fname)
load(object_fname_pth)

#for file names
dbCO="top10"

#how many dbs identify the pair (for predicted miRNA-mRNA pairs)
ndb=1

readme_fname=file.path(resdir,"readme")
writeLines(c(
  paste("mRNA","FDR",FDR_CO, sep="\t"),
  paste("mRNA","log2FC",lfc_CO, sep="\t") ,
  paste("miRNA","FDR",FDR_COmir, sep="\t"),
  paste("miRNA","log2FC",lfc_COmir, sep="\t") ,
  paste("predicted miRNA - mRNA targets: top",dbCO,"db hits (ranked by p value)", sep=" ") ,
  paste("number of dbs which predict miRNA - mRNA target:",ndb, sep=" ") ),
   readme_fname)



# rbind summaries to obtain one df for targets in each class (2 x predicted and validated)

multimir_summary_validated=list(all_genes_array_targeted_miRNA_valid_1@summary,
  all_genes_array_targeted_miRNA_valid_2@summary,
  all_genes_array_targeted_miRNA_valid_3@summary)


all_genes_array_targeted_miRNA_v=tibble(data.table::rbindlist(multimir_summary_validated, fill=TRUE, idcol=NULL))

# > nrow(all_genes_array_targeted_miRNA_v)
# [1] 121926

# two sets of predicted targets, top10

multimir_summary_pred_cons=list(all_genes_array_targeted_miRNA_pred_1@summary,
  all_genes_array_targeted_miRNA_pred_2@summary,
  all_genes_array_targeted_miRNA_pred_3@summary)



all_genes_array_targeted_miRNA_p_cons=tibble(data.table::rbindlist(multimir_summary_pred_cons, fill=TRUE, idcol=NULL))

# > nrow(all_genes_array_targeted_miRNA_p_cons)
# [1] 170580


multimir_summary_pred_all=list(a10_all_genes_array_targeted_miRNA_pred_1@summary,
  a10_all_genes_array_targeted_miRNA_pred_2@summary,
  a10_all_genes_array_targeted_miRNA_pred_3@summary,
  a10_all_genes_array_targeted_miRNA_pred_4@summary,
  a10_all_genes_array_targeted_miRNA_pred_5@summary)

all_genes_array_targeted_miRNA_p_all=tibble(data.table::rbindlist(multimir_summary_pred_all, fill=TRUE, idcol=NULL))

# > nrow(all_genes_array_targeted_miRNA_p_all)
# [1] 346787



# two sets of predicted targets, top20

# multimir_summary_pred_cons=list(c20_all_genes_array_targeted_miRNA_pred_1@summary,
#   c20_all_genes_array_targeted_miRNA_pred_2@summary,
#   c20_all_genes_array_targeted_miRNA_pred_3@summary)



# all_genes_array_targeted_miRNA_p_cons=tibble(data.table::rbindlist(multimir_summary_pred_cons, fill=TRUE, idcol=NULL))

# > nrow(all_genes_array_targeted_miRNA_p_cons)
# [1] 170580


# multimir_summary_pred_all=list(a20_all_genes_array_targeted_miRNA_pred_1@summary,
#   a20_all_genes_array_targeted_miRNA_pred_2@summary,
#   a20_all_genes_array_targeted_miRNA_pred_3@summary,
#   a20_all_genes_array_targeted_miRNA_pred_4@summary,
#   a20_all_genes_array_targeted_miRNA_pred_5@summary)

# all_genes_array_targeted_miRNA_p_all=tibble(data.table::rbindlist(multimir_summary_pred_all, fill=TRUE, idcol=NULL))



## ---- get_targets_subgroups


# the targets are retrieved from the saved multimir objects


# #mRNA
# mRNA_up=results_mRNA_sig[results_mRNA_sig$logFC>lfc_CO,]
# mRNA_dn=results_mRNA_sig[results_mRNA_sig$logFC<(-lfc_CO),]


# #miRNA
# mir_up=mir_sig[mir_sig$logFC>lfc_COmir,]
# mir_dn=mir_sig[mir_sig$logFC<(-lfc_COmir),]




# mir UP genes DN
mirnas=rownames(mir_up)
genes=mRNA_dn$SYMBOL

#top10 or top20

miRup_GENdn_pred_cons=select(all_genes_array_targeted_miRNA_p_cons, 
  mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)


miRup_GENdn_pred_all=select(all_genes_array_targeted_miRNA_p_all,
   mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)


miRup_GENdn_valid=select(all_genes_array_targeted_miRNA_v,
 mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,validated.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)


# mir UP genes UP
mirnas=rownames(mir_up)
genes=mRNA_up$SYMBOL

miRup_GENup_pred_cons=select(all_genes_array_targeted_miRNA_p_cons, 
  mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)

miRup_GENup_pred_all=select(all_genes_array_targeted_miRNA_p_all,
   mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)


miRup_GENup_valid=select(all_genes_array_targeted_miRNA_v,
 mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,validated.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)


# mir DN genes DN
mirnas=rownames(mir_dn)
genes=mRNA_dn$SYMBOL


miRdn_GENdn_pred_cons=select(all_genes_array_targeted_miRNA_p_cons, 
  mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)

miRdn_GENdn_pred_all=select(all_genes_array_targeted_miRNA_p_all,
   mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)


miRdn_GENdn_valid=select(all_genes_array_targeted_miRNA_v,
 mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,validated.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)

# mir DN genes UP
mirnas=rownames(mir_dn)
genes=mRNA_up$SYMBOL

miRdn_GENup_pred_cons=select(all_genes_array_targeted_miRNA_p_cons, 
  mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)

miRdn_GENup_pred_all=select(all_genes_array_targeted_miRNA_p_all,
   mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,predicted.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)%>%
  filter(predicted.sum>=ndb)


miRdn_GENup_valid=select(all_genes_array_targeted_miRNA_v,
 mature_mirna_acc,mature_mirna_id,target_symbol,target_entrez,target_ensembl,validated.sum) %>%
  filter(mature_mirna_id%in%mirnas) %>%
  filter(target_symbol%in%genes)






## ---- summary_retrieved_numbers
#summary of identified targets
# be mindful of top10 or top20

targets_sum=data.frame(Validated=c(nrow(miRup_GENdn_valid),nrow(miRup_GENup_valid),nrow(miRdn_GENdn_valid),nrow(miRdn_GENup_valid) ),
   Predicted_conserved=c(nrow(miRup_GENdn_pred_cons),nrow(miRup_GENup_pred_cons),nrow(miRdn_GENdn_pred_cons),nrow(miRdn_GENup_pred_cons) ),
   Predicted_all=c(nrow(miRup_GENdn_pred_all),nrow(miRup_GENup_pred_all),nrow(miRdn_GENdn_pred_all),nrow(miRdn_GENup_pred_all) ) )
rownames(targets_sum)=c("miRup_GENdn","miRup_GENup","miRdn_GENdn","miRdn_GENup")


#logFC=0, FDR.mRNA=0.05, FDR.miRNA=0.01
# top10 conserved sites
#             Validated Predicted
# miRup_GENdn       116       118
# miRup_GENup       170       205
# miRdn_GENdn       268       310
# miRdn_GENup       411       454


# > lfc_CO
# [1] 0
# > lfc_COmir
# [1] 0
# > FDR_CO
# [1] 0.1
# > FDR_COmir
# [1] 0.01
#             Validated Predicted_top10_cons Predicted_top20_all
# miRup_GENdn       403                  481                2324
# miRup_GENup       507                  671                2599
# miRdn_GENdn       940                 1076                4269
# miRdn_GENup      1450                 1686                5100



#table for the report
caption.tab=c("Number of mRNA - miRNA interactions in each expression category. 
	For example \\texttt{miRup-GENdn} are targets of up-regulated miRNAs (up-072 data set) 
	which are amongst genes downregulated in the mRNA expression profiling (GSE12049 data set).")


print(xtable(as.data.frame(targets_sum) , align=c("p{4cm}","c","c","c"), 
display=c("s","d","d","d"),
caption=caption.tab),
include.rownames=TRUE,floating=TRUE,
table.placement="H", caption.placement = "bottom",
size="\\scriptsize")





## ---- TEST_list_overlap_notes_only



#background sets
# all records related to miRNAs in the experiment and genes on the array (after filtering) (in predicted and validated dbs)

# > length(unique(all_genes_array_targeted_miRNA_p_all$target_symbol))
# [1] 13344
# > 
# > length(unique(all_genes_array_targeted_miRNA_p_cons$target_symbol))
# [1] 12128
# > 
# > length(unique(all_genes_array_targeted_miRNA_v$target_symbol))
# [1] 12236


## ---- function_test_overlap

# this section is to test whether there is a sginificant overlap between (i) the list of DE genes and 
# (ii) the list of miRNA targets (predicted and validated) for each queried miRNA


# mir_to_test - list of miRNAs to test (from df of significant DE results)

# mRNA_all - all mRNAs on the array
# mRNA_all=results_mRNA$SYMBOL

# mRNA_changed - mRNA_dn$SYMBOL; list of up / dn mRNAs to be tested
# mRNA_changed=length(mRNA_dn$SYMBOL)

# targets.v, targets.p - data frames with summary from multiMir query results
# this gives changed genes targetd by SELECTED miRNA GROUP in this experiment
# targetsV_changed=as.data.frame(miRup_GENdn_valid@summary)
# targetsP_changed=as.data.frame(miRup_GENdn_pred@summary)
# data frame
#head(targetsV_changed)
#   mature_mirna_acc mature_mirna_id target_symbol target_entrez     target_ensembl
# 1     MIMAT0000154 mmu-miR-142a-5p          Cd28         12487 ENSMUSG00000026012
# 2     MIMAT0000162  mmu-miR-152-3p        Camk2a         12322 ENSMUSG00000024617
# 3     MIMAT0000527   mmu-miR-16-5p         Vegfa         22339 ENSMUSG00000023951
# 4     MIMAT0000661  mmu-miR-214-3p         Cd59a         12509 ENSMUSG00000032679
# 5     MIMAT0000665  mmu-miR-223-3p        Cd99l2        171486 ENSMUSG00000035776
# 6     MIMAT0000665  mmu-miR-223-3p         Epdr1        105298 ENSMUSG00000002808
#   mirtarbase tarbase validated.sum all.sum
# 1          2       1             2       2
# 2          1       1             2       2


# this gives all genes targetd by SELECTED miRNA in this experiment
# targetsV_all=all_genes_array_targeted_miRNA_valid
# targetsP_all=all_genes_array_targeted_miRNA_pred
# data frame (or tibble)
#   mature_mirna_acc mature_mirna_id target_symbol target_entrez target_ensembl    
#   <chr>            <chr>           <chr>         <chr>         <chr>             
# 1 MIMAT0000123     mmu-miR-1a-3p   Cdc42         12540         ENSMUSG00000006699



#this function returns table with results of overlap of lists of targets of selected (DE) miRNA targets and DE mRNAs
# input are targets for predicted conserved, predicted all, validated etc.


test_overlap <- function(mir_to_test,mRNA_changed,mRNA_all,targets_changed,targets_all){

  DE_genes=as.numeric(length(unique(mRNA_changed))) # count of DE genes

  all_genes_array=as.numeric(length(unique(mRNA_all))) # all genes on the array (background set)

  targets_changed=as.data.frame(targets_changed)
  targets_all=as.data.frame(targets_all)


  test_res=NULL
  test_res=data.frame(miRNA=as.character(),pval=as.numeric(),
    n_mRNA_DE=as.character(),n_targets_DE=as.character(),
    n_targets=as.character(),n_mRNA=as.character())

  for (i in mir_to_test){


    # all genes in miRNA-target db associated with this miRNA
    genes_array_miRNA.i=targets_all[targets_all$mature_mirna_id==i,]
    genes_targ_miR.i=as.numeric(length(unique(genes_array_miRNA.i$target_symbol)))

    #DE genes in miRNA-target db associated with this miRNA
    targets.i=targets_changed[targets_changed$mature_mirna_id==i,]
    DE_targets=as.numeric(length(unique(targets.i$target_symbol)))

    # test for list overlap
    overlap=DE_targets
    ALLtargets.miRNA.i=genes_targ_miR.i

    #print(is(DE_genes))
    #print(is(all_genes_array))
    #print(is(ALLtargets.miRNA.i))
    #print(is(overlap))

    p.hyper=phyper(overlap-1 , DE_genes , all_genes_array-DE_genes , ALLtargets.miRNA.i ,lower.tail=FALSE, log.p=FALSE)

    # p.fisher=fisher.test(matrix(c( 
    #   overlap , DE_genes-overlap , ALLtargets.miRNA.i-overlap , all_genes_array-DE_genes-ALLtargets.miRNA.i+overlap),
    #   2,2),  alternative='greater')$p.value
    
    miR_line=c(i, p.hyper,DE_targets,genes_targ_miR.i,DE_genes,all_genes_array)

    test_res=rbind(test_res,miR_line)

}

  colnames(test_res)=c("miRNA","p.hyper","n_DE_targets","n_targets","n_DE_genes","all_genes")
  #test_res$p.hyper=formatC(as.numeric(test_res$p.hyper), format="fg")
  test_res$p.hyper=formatC(as.numeric(test_res$p.hyper), format="g")


  return(test_res)
}

#wrapper for the above which also formats the results
# input are data frames: in pairs all and changed
# output is a list of dfs: with results of overlap (formatted) and list of genes

get_overlaps<-function(targets_all,targets_changed,name,mir_to_test,mRNA_all,mRNA_changed){
  
  list_overlap_res=NULL
  list_overlap_res=list()

  overlap_res=test_overlap(mir_to_test=mir_to_test,mRNA_all=mRNA_all,
  mRNA_changed=mRNA_changed,targets_changed=targets_changed,targets_all=targets_all)

  gene_totals=overlap_res[,c(1,5,6)]

  colnames(overlap_res)[2:ncol(overlap_res)]=paste(colnames(overlap_res)[2:ncol(overlap_res)],name, sep=".") #formatting for cbind later
  overlap_results_v=overlap_res[,-c(5,6)]

  list_overlap_res[[1]]=overlap_results_v
  list_overlap_res[[2]]=gene_totals

  return(list_overlap_res)
}


save_df_overlaps<-function(list_overlap_results,fpth){

  ## final df
  overlap_results_df = rlist::list.cbind(list_overlap_results)
  rownames(overlap_results_df) = overlap_results_df$miRNA
  overlap_results_df = overlap_results_df %>% select(-contains("miRNA"))
  overlap_results_df$miRNA = rownames(overlap_results_df)
  overlap_results_df = overlap_results_df %>% select(miRNA, everything())

  write.table(overlap_results_df,fpth,sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  return(overlap_results_df)
}


## ---- TEST_overlap_not_in_report

# example usage

mRNA_all=unique(results_mRNA$SYMBOL)
mRNA_changed=unique(mRNA_up$SYMBOL)
mir_to_test=mir_dn$miR
targets_changed=miRdn_GENup_valid
targets_all=all_genes_array_targeted_miRNA_v


results_tst=test_overlap(mir_to_test=mir_to_test,mRNA_all=mRNA_all,
  mRNA_changed=mRNA_changed,targets_changed=targets_changed,targets_all=targets_all)


#########
# comparisons

## ---- TEST_list_overlaps



mi_mRNA_int.dirname=paste("miRNA-mRNA.list-overlap",dbCO,sep=".")
res_miRNA_mRNA_overlap=file.path(resdir,mi_mRNA_int.dirname)
dir.create(res_miRNA_mRNA_overlap)


##############################
# miRdn_GENup

fname=paste("ListOverlap.miRNA_DN.mRNA_UP",dbCO,"tsv",sep=".")
fpth=file.path(res_miRNA_mRNA_overlap,fname)

mRNA_all=unique(results_mRNA$SYMBOL)

mir_to_test=mir_dn$miR
mRNA_changed=unique(mRNA_up$SYMBOL)

list_overlap_results=NULL
list_overlap_results=list()

rv=get_overlaps(targets_all=all_genes_array_targeted_miRNA_v,targets_changed=miRdn_GENup_valid,name="validated",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[1]]=rv[[1]]

rpc=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_cons,targets_changed=miRdn_GENup_pred_cons,name="predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[2]]=rpc[[1]]

rpa=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_all,targets_changed=miRdn_GENup_pred_all,name="predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[3]]=rpa[[1]]


#COMBINED predicted + validated

#conserved

targets_all=rbind(all_genes_array_targeted_miRNA_p_cons[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRdn_GENup_pred_cons[,1:5],
  miRdn_GENup_valid[,1:5])


rvpc=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[4]]=rvpc[[1]]

#all

targets_all=rbind(all_genes_array_targeted_miRNA_p_all[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRdn_GENup_pred_all[,1:5],
  miRdn_GENup_valid[,1:5])

rvpa=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[5]]=rvpa[[1]]

#totals to add to the final df
list_overlap_results[[6]]=rvpa[[2]]

overlap_results_df=save_df_overlaps(list_overlap_results=list_overlap_results,fpth=fpth)

sig.miRdn_GENup.predCons=overlap_results_df[overlap_results_df$p.hyper.predCons<0.05,]
sig.miRdn_GENup.predAll=overlap_results_df[overlap_results_df$p.hyper.predAll<0.05,]
sig.miRdn_GENup.predCons.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predCons<0.05,]
sig.miRdn_GENup.predAll.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predAll<0.05,]



##############################
# miRup_GENdn

fname=paste("ListOverlap.miRNA_UP.mRNA_DN",dbCO,"tsv",sep=".")
fpth=file.path(res_miRNA_mRNA_overlap,fname)


mir_to_test=mir_up$miR
mRNA_changed=unique(mRNA_dn$SYMBOL)

list_overlap_results=NULL
list_overlap_results=list()

rv=get_overlaps(targets_all=all_genes_array_targeted_miRNA_v,targets_changed=miRup_GENdn_valid,name="validated",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[1]]=rv[[1]]

rpc=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_cons,targets_changed=miRup_GENdn_pred_cons,name="predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[2]]=rpc[[1]]

rpa=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_all,targets_changed=miRup_GENdn_pred_all,name="predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[3]]=rpa[[1]]


#COMBINED predicted + validated

#conserved

targets_all=rbind(all_genes_array_targeted_miRNA_p_cons[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRup_GENdn_pred_cons[,1:5],
  miRup_GENdn_valid[,1:5])


rvpc=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[4]]=rvpc[[1]]

#all

targets_all=rbind(all_genes_array_targeted_miRNA_p_all[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRup_GENdn_pred_all[,1:5],
  miRup_GENdn_valid[,1:5])

rvpa=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[5]]=rvpa[[1]]

#totals to add to the final df
list_overlap_results[[6]]=rvpa[[2]]

overlap_results_df=save_df_overlaps(list_overlap_results=list_overlap_results,fpth=fpth)

sig.miRup_GENdn.predCons=overlap_results_df[overlap_results_df$p.hyper.predCons<0.05,]
sig.miRup_GENdn.predAll=overlap_results_df[overlap_results_df$p.hyper.predAll<0.05,]
sig.miRup_GENdn.predCons.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predCons<0.05,]
sig.miRup_GENdn.predAll.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predAll<0.05,]



##############################
# miRdn_GENdn

fname=paste("ListOverlap.miRNA_DN.mRNA_DN",dbCO,"tsv",sep=".")
fpth=file.path(res_miRNA_mRNA_overlap,fname)


mir_to_test=mir_dn$miR
mRNA_changed=unique(mRNA_dn$SYMBOL)

list_overlap_results=NULL
list_overlap_results=list()

rv=get_overlaps(targets_all=all_genes_array_targeted_miRNA_v,targets_changed=miRdn_GENdn_valid,name="validated",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[1]]=rv[[1]]

rpc=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_cons,targets_changed=miRdn_GENdn_pred_cons,name="predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[2]]=rpc[[1]]

rpa=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_all,targets_changed=miRdn_GENdn_pred_all,name="predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[3]]=rpa[[1]]


#COMBINED predicted + validated

#conserved

targets_all=rbind(all_genes_array_targeted_miRNA_p_cons[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRdn_GENdn_pred_cons[,1:5],
  miRdn_GENdn_valid[,1:5])


rvpc=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[4]]=rvpc[[1]]

#all

targets_all=rbind(all_genes_array_targeted_miRNA_p_all[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRdn_GENdn_pred_all[,1:5],
  miRdn_GENdn_valid[,1:5])

rvpa=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[5]]=rvpa[[1]]

#totals to add to the final df
list_overlap_results[[6]]=rvpa[[2]]

overlap_results_df=save_df_overlaps(list_overlap_results=list_overlap_results,fpth=fpth)

sig.miRdn_GENdn.predCons=overlap_results_df[overlap_results_df$p.hyper.predCons<0.05,]
sig.miRdn_GENdn.predAll=overlap_results_df[overlap_results_df$p.hyper.predAll<0.05,]
sig.miRdn_GENdn.predCons.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predCons<0.05,]
sig.miRdn_GENdn.predAll.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predAll<0.05,]



##############################
# miRup_GENup

fname=paste("ListOverlap.miRNA_UP.mRNA_UP",dbCO,"tsv",sep=".")
fpth=file.path(res_miRNA_mRNA_overlap,fname)


mir_to_test=mir_up$miR
mRNA_changed=unique(mRNA_up$SYMBOL)

list_overlap_results=NULL
list_overlap_results=list()

rv=get_overlaps(targets_all=all_genes_array_targeted_miRNA_v,targets_changed=miRup_GENup_valid,name="validated",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[1]]=rv[[1]]

rpc=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_cons,targets_changed=miRup_GENup_pred_cons,name="predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[2]]=rpc[[1]]

rpa=get_overlaps(targets_all=all_genes_array_targeted_miRNA_p_all,targets_changed=miRup_GENup_pred_all,name="predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[3]]=rpa[[1]]


#COMBINED predicted + validated

#conserved

targets_all=rbind(all_genes_array_targeted_miRNA_p_cons[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRup_GENup_pred_cons[,1:5],
  miRup_GENup_valid[,1:5])


rvpc=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predCons",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[4]]=rvpc[[1]]

#all

targets_all=rbind(all_genes_array_targeted_miRNA_p_all[,1:5], 
  all_genes_array_targeted_miRNA_v[,1:5])

targets_changed=rbind(miRup_GENup_pred_all[,1:5],
  miRup_GENup_valid[,1:5])

rvpa=get_overlaps(targets_all=targets_all,targets_changed=targets_changed,name="combined.valid.predAll",
  mir_to_test=mir_to_test,mRNA_all=mRNA_all,mRNA_changed=mRNA_changed)

list_overlap_results[[5]]=rvpa[[1]]

#totals to add to the final df
list_overlap_results[[6]]=rvpa[[2]]

overlap_results_df=save_df_overlaps(list_overlap_results=list_overlap_results,fpth=fpth)

sig.miRup_GENup.predCons=overlap_results_df[overlap_results_df$p.hyper.predCons<0.05,]
sig.miRup_GENup.predAll=overlap_results_df[overlap_results_df$p.hyper.predAll<0.05,]
sig.miRup_GENup.predCons.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predCons<0.05,]
sig.miRup_GENup.predAll.comb=overlap_results_df[overlap_results_df$p.hyper.combined.valid.predAll<0.05,]




## ---- summary_sig_miRs

summary_sig_miRs=data.frame(Comparison=c("miRNA UP, mRNA DOWN","miRNA DOWN, mRNA DOWN","miRNA UP, mRNA UP","miRNA DOWN, mRNA UP"),
	predicted_conserved=c(nrow(sig.miRup_GENdn.predCons),nrow(sig.miRdn_GENdn.predCons),nrow(sig.miRup_GENup.predCons),nrow(sig.miRdn_GENup.predCons)),
  combined_predicted_conserved_validated=c(nrow(sig.miRup_GENdn.predCons.comb),nrow(sig.miRdn_GENdn.predCons.comb),nrow(sig.miRup_GENup.predCons.comb),nrow(sig.miRdn_GENup.predCons.comb)),
  predicted_all=c(nrow(sig.miRup_GENdn.predAll),nrow(sig.miRdn_GENdn.predAll),nrow(sig.miRup_GENup.predAll),nrow(sig.miRdn_GENup.predAll)),
  combined_predicted_all_validated=c(nrow(sig.miRup_GENdn.predAll.comb),nrow(sig.miRdn_GENdn.predAll.comb),nrow(sig.miRup_GENup.predAll.comb),nrow(sig.miRdn_GENup.predAll.comb))
)


colnames(summary_sig_miRs)=c("Comparison","Targets Predicted Conserved","Targets Combined Validated + Predicted Conserved","Targets Predicted All","Targets Combined Validated + Predicted All")




#table for the report
caption.tab=c("Number of significant overlaps of differentially expressed mRNA genes with targets of differentially expressed miRNAs. 
	For example \\texttt{miRNA UP, mRNA DOWN} are targets of up-regulated miRNAs (up-072 data set) 
	which are amongst genes downregulated in the mRNA expression profiling (GSE12049 data set).")


print(xtable(as.data.frame(summary_sig_miRs) , align=c("p{0.1cm}","c","p{2.8cm}","p{2.8cm}","p{2.8cm}","p{2.8cm}"), 
display=c("s","s","d","d","d","d"),
caption=caption.tab),
include.rownames=FALSE,floating=TRUE,
table.placement="H", caption.placement = "bottom",size="\\scriptsize")


#######################################
#######################################
#######################################
#######################################
#######################################
#######################################
#######################################


## ---- sig_miRs_details

#xtables from all expression status groups

sig.miRdn_GENup=rbind(sig.miRdn_GENup.predCons,sig.miRdn_GENup.predAll,sig.miRdn_GENup.predCons.comb,sig.miRdn_GENup.predAll.comb)%>% distinct()

table_rep=sig.miRdn_GENup[,c(1,2,5,8,11,14)]
colnames(table_rep)=c("miRNA","p.hyper validated","p.hyper predicted conserved","p.hyper predicted all","p.hyper predicted conserved + validated","p.hyper predicted all + validated")

caption.tab=c("Down-regulated miRNAs with significant overlaps of their targets 
  predicted (conserved and all sites) and validated with up-regulated mRNA genes (GSE12049 data set).")

nrows=nrow(table_rep)

if (nrows >0){

  if (nrows >1){ 
    rws=seq(1, (nrow(table_rep)-1), by = 2)
    col=rep("\\rowcolor[gray]{0.95}", length(rws))
    }else{
      rws=1
      col=rep("\\rowcolor[gray]{0.95}", length(rws))
  }

  print(xtable(as.data.frame(table_rep) , align=c("p{0.1cm}","c","p{2cm}","p{2cm}","p{2cm}","p{2cm}","p{2cm}"), 
  display=c("s","s","e","e","e","e","e"),  #digits=c(0,0,3,3,3),
  caption=caption.tab),
  include.rownames=FALSE,floating=TRUE,
  table.placement="H", caption.placement = "bottom",
  add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")

}

sig.miRup_GENdn=rbind(sig.miRup_GENdn.predCons,sig.miRup_GENdn.predAll,sig.miRup_GENdn.predCons.comb,sig.miRup_GENdn.predAll.comb)%>% distinct()

table_rep=sig.miRup_GENdn[,c(1,2,5,8,11,14)]
colnames(table_rep)=c("miRNA","p.hyper validated","p.hyper predicted conserved","p.hyper predicted all","p.hyper predicted conserved + validated","p.hyper predicted all + validated")

caption.tab=c("Up-regulated miRNAs with significant overlaps of their targets 
  predicted (conserved and all sites) and validated with down-regulated mRNA genes (GSE12049 data set).")

nrows=nrow(table_rep)

if (nrows >0){

  if (nrows >1){ 
    rws=seq(1, (nrow(table_rep)-1), by = 2)
    col=rep("\\rowcolor[gray]{0.95}", length(rws))
    }else{
      rws=1
      col=rep("\\rowcolor[gray]{0.95}", length(rws))
  }

  print(xtable(as.data.frame(table_rep) , align=c("p{0.1cm}","c","p{2cm}","p{2cm}","p{2cm}","p{2cm}","p{2cm}"), 
  display=c("s","s","e","e","e","e","e"),  #digits=c(0,0,3,3,3),
  caption=caption.tab),
  include.rownames=FALSE,floating=TRUE,
  table.placement="H", caption.placement = "bottom",
  add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")

}

sig.miRup_GENup=rbind(sig.miRup_GENup.predCons,sig.miRup_GENup.predAll,sig.miRup_GENup.predCons.comb,sig.miRup_GENup.predAll.comb)%>% distinct()

table_rep=sig.miRup_GENup[,c(1,2,5,8,11,14)]
colnames(table_rep)=c("miRNA","p.hyper validated","p.hyper predicted conserved","p.hyper predicted all","p.hyper predicted conserved + validated","p.hyper predicted all + validated")

caption.tab=c("Up-regulated miRNAs with significant overlaps of their targets 
  predicted (conserved and all sites) and validated with up-regulated mRNA genes (GSE12049 data set).")

nrows=nrow(table_rep)

if (nrows >0){

  if (nrows >1){ 
    rws=seq(1, (nrow(table_rep)-1), by = 2)
    col=rep("\\rowcolor[gray]{0.95}", length(rws))
    }else{
      rws=1
      col=rep("\\rowcolor[gray]{0.95}", length(rws))
  }

  print(xtable(as.data.frame(table_rep) , align=c("p{0.1cm}","c","p{2cm}","p{2cm}","p{2cm}","p{2cm}","p{2cm}"), 
  display=c("s","s","e","e","e","e","e"),  #digits=c(0,0,3,3,3),
  caption=caption.tab),
  include.rownames=FALSE,floating=TRUE,
  table.placement="H", caption.placement = "bottom",
  add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")

}

sig.miRdn_GENdn=rbind(sig.miRdn_GENdn.predCons,sig.miRdn_GENdn.predAll,sig.miRdn_GENdn.predCons.comb,sig.miRdn_GENdn.predAll.comb)%>% distinct()

table_rep=sig.miRdn_GENdn[,c(1,2,5,8,11,14)]
colnames(table_rep)=c("miRNA","p.hyper validated","p.hyper predicted conserved","p.hyper predicted all","p.hyper predicted conserved + validated","p.hyper predicted all + validated")

caption.tab=c("Down-regulated miRNAs with significant overlaps of their targets 
  predicted (conserved and all sites) and validated with down-regulated mRNA genes (GSE12049 data set).")

nrows=nrow(table_rep)

if (nrows >0){

  if (nrows >1){ 
    rws=seq(1, (nrow(table_rep)-1), by = 2)
    col=rep("\\rowcolor[gray]{0.95}", length(rws))
    }else{
      rws=1
      col=rep("\\rowcolor[gray]{0.95}", length(rws))
  }

  print(xtable(as.data.frame(table_rep) , align=c("p{0.1cm}","c","p{2cm}","p{2cm}","p{2cm}","p{2cm}","p{2cm}"), 
  display=c("s","s","e","e","e","e","e"),  #digits=c(0,0,3,3,3),
  caption=caption.tab),
  include.rownames=FALSE,floating=TRUE,
  table.placement="H", caption.placement = "bottom",
  add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")

}


## ---- sel_miRNA_prep

results_mRNA_short=results_mRNA[,c(1,2,3,4,5,6,8,9)]

# for these visualisations select best probe for each gene; best i.e. significant and largest abs(logFC)

results_mRNA_short_sig=results_mRNA_short%>%
  filter(adj.P.Val<FDR_CO)%>%
  filter(abs(logFC)>lfc_CO)%>%
  group_by(SYMBOL)%>%
  slice(which.max(abs(logFC)))

# # A tibble: 675 x 8
# # Groups:   SYMBOL [675]
#    PROBEID      SYMBOL        GENENAME                   ENSEMBL             logFC AveExpr    P.Value adj.P.Val
#    <chr>        <chr>         <chr>                      <chr>               <dbl>   <dbl>      <dbl>     <dbl>
#  1 1457779_at   1110046J04Rik RIKEN cDNA 1110046J04 gene ENSMUSG00000085457  1.08     4.59 0.00117      0.0707 
#  2 1431084_x_at 1810014B01Rik RIKEN cDNA 1810014B01 gene ENSMUSG00000097412 -0.892    5.68 0.00111      0.0696 
#  3 1432420_a_at 2310002L09Rik RIKEN cDNA 2310002L09 gene ENSMUSG00000028396 -1.27     8.52 0.000258     0.0445 
#  4 1431572_at   2310015K22Rik RIKEN cDNA 2310015K22 gene NA                 -1.38     5.23 0.000501     0.0501 
#  5 1443575_at   2310040G24Rik RIKEN cDNA 2310040G24 gene ENSMUSG00000101655 -1.46     8.74 0.000535     0.0511 
#  6 1440008_at   2310043L19Rik RIKEN cDNA 2310043L19 gene ENSMUSG00000101746  1.08     4.56 0.000275     0.0445 
#  7 1458587_at   2310047D07Rik RIKEN cDNA 2310047D07 gene NA                 -1.09     6.61 0.00111      0.0696 
#  8 1453657_at   2310065F04Rik RIKEN cDNA 2310065F04 gene ENSMUSG00000087410 -2.69     7.00 0.00000434   0.00875
#  9 1447936_at   2410006H16Rik RIKEN cDNA 2410006H16 gene ENSMUSG00000086841  1.18     8.17 0.00108      0.0694 
# 10 1438429_at   2610319H10Rik RIKEN cDNA 2610319H10 gene NA                  1.18     7.54 0.00312      0.0986 
# # … with 665 more rows

#annot for hm
colour_palette_hm=colorRamp2(c(-2, 0, 2), c("dodgerblue4", "ivory", "chocolate3")) #z-scored

annot_col=samples
annot_col=as.data.frame(annot_col[,-1])
colnames(annot_col)="Strain"

ann_colors=list(
  Strain=c(Mut="darkorange2", WT="dodgerblue4"))


## ---- example_sig_mir_report

# this section is to prep data for tables and plotting for a selected example miRNA
i="mmu-miR-1a-3p"

caption_levels.x=paste("Median normalised log expression of\n",i,sep=" ")

# pred.cons=miRdn_GENup_pred_cons
# pred.all=miRdn_GENup_pred_all
# valid=miRdn_GENup_valid

#sometimes the same target_symbol - miR combination appears more than once
# > pred.cons[pred.cons$target_symbol=="Matr3",]
# # A tibble: 16 x 6
#    mature_mirna_acc mature_mirna_id target_symbol target_entrez target_ensembl     predicted.sum
#    <chr>            <chr>           <chr>         <chr>         <chr>                      <dbl>
#  1 MIMAT0000123     mmu-miR-1a-3p   Matr3         "17184"       ENSMUSG00000037236             4
#  2 MIMAT0000239     mmu-miR-206-3p  Matr3         "17184"       ENSMUSG00000037236             4
#  3 MIMAT0000123     mmu-miR-1a-3p   Matr3         ""            ENSMUSG00000037236             2
#remove the additional occurences to keep only one

pred.cons=miRdn_GENup_pred_cons%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

pred.all=miRdn_GENup_pred_all%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

valid=miRdn_GENup_valid


pred.cons.i=pred.cons[pred.cons$mature_mirna_id==i,]
pred.all.i=pred.all[pred.all$mature_mirna_id==i,]
valid.i=valid[valid$mature_mirna_id==i,]

pred.cons.i_ann=left_join(pred.cons.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
  filter(adj.P.Val<FDR_CO) %>%
  filter(abs(logFC)>lfc_CO)

pred.all.i_ann=left_join(pred.all.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
  filter(adj.P.Val<FDR_CO) %>%
  filter(abs(logFC)>lfc_CO)

valid.i_ann=left_join(valid.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
  filter(adj.P.Val<FDR_CO) %>%
  filter(abs(logFC)>lfc_CO)

# > pred.cons.i_ann[,c(2,3,6,7,8,10,13)]
# # A tibble: 39 x 7
#    mature_mirna_id target_symbol predicted.sum PROBEID  GENENAME       logFC adj.P.Val
#    <chr>           <chr>                 <dbl> <chr>    <chr>          <dbl>     <dbl>
#  1 mmu-miR-1a-3p   Cnn3                      5 1436759… calponin 3, a… 0.926    0.0870
#  2 mmu-miR-1a-3p   Foxp1                     4 1435222… forkhead box … 0.909    0.0513
#  3 mmu-miR-1a-3p   Matr3                     4 1438368… matrin 3       0.895    0.0845
#  4 mmu-miR-1a-3p   Ncl                       4 1456528… nucleolin      0.672    0.0906
#  5 mmu-miR-1a-3p   Hnrnpu                    3 1423050… heterogeneous… 0.631    0.0920
#  6 mmu-miR-1a-3p   Sptbn1                    3 1444089… spectrin beta… 1.04     0.0501
#  7 mmu-miR-1a-3p   Ywhab                     3 1420880… tyrosine 3-mo… 0.629    0.0886
#  8 mmu-miR-1a-3p   Abca1                     2 1421840… ATP-binding c… 1.20     0.0859
#  9 mmu-miR-1a-3p   Eid1                      2 1448405… EP300 interac… 0.846    0.0450
# 10 mmu-miR-1a-3p   Eid1                      2 1416614… EP300 interac… 0.762    0.0853


##### visualisation

## exp miR i
expr_miRNA.i=as.numeric(exp_mir[rownames(exp_mir)==i,])

## median of the expr for plotting
miRNA_median.wt=median(expr_miRNA.i[6:10])
miRNA_median.mut=median(expr_miRNA.i[1:5])


## top 3 targets for vis (predicted conserved) -  ranked by log FC
targets.ann=pred.cons.i_ann

targets.ann_sel=targets.ann[order(abs(targets.ann$logFC), decreasing=TRUE),][1:3,]
targets.ann_ord_names=targets.ann_sel$target_symbol
targets.ann_ord_probes=targets.ann_sel$PROBEID
names(targets.ann_ord_probes)=targets.ann_ord_names


#loop_for_plots
list_of_plots=NULL
list_of_captions=NULL


list_to_plot=targets.ann_ord_probes

n_probs=length(list_to_plot)

for (jj in (1:n_probs)){

  j=list_to_plot[jj]
  name.j=names(list_to_plot[jj])

  caption_levels.y=paste("Normalised log expression of\n",name.j,sep=" ")

  mRNA_exprs.j=exprs(data_final[rownames(data_final)==j,])

  df_plot.j=as.data.frame(t(mRNA_exprs.j))
  colnames(df_plot.j)="mRNA"

  df_plot.j$strain=c(rep("WT",3),rep("Mut",3))

  df_plot.j$miRNA_median=c(rep(miRNA_median.wt,3),rep(miRNA_median.mut,3))

  plot.j_save=ggplot(df_plot.j, aes(x=miRNA_median, y=mRNA,colour=strain)) + geom_boxplot(width=0.5) + geom_point() +
  scale_color_manual(values = c("darkorange2", "dodgerblue4")) + theme_pca +
  xlab(caption_levels.x)+ ylab(caption_levels.y)

  plot.j=ggplot(df_plot.j, aes(x=miRNA_median, y=mRNA,colour=strain)) + geom_boxplot(width=0.5) + geom_point() +
  scale_color_manual(values = c("darkorange2", "dodgerblue4")) + theme_pca +
  xlab(caption_levels.x)+ ylab(caption_levels.y) + fontsize2

  list_of_plots[[jj]]=plot.j

  caption.j=paste("Expression of",name.j,"vs. median expression of",i,"in WT and Mut.")

  list_of_captions[[jj]]=caption.j

}


## ---- heatmaps_sel_mir

pred_GENup=miRdn_GENup_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
pred_GENdn=miRdn_GENdn_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))

pred_GENup.i=pred_GENup[pred_GENup$mature_mirna_id==i,]
pred_GENdn.i=pred_GENdn[pred_GENdn$mature_mirna_id==i,]

pred_all=rbind(pred_GENup.i,pred_GENdn.i)

# for this heatmap select best probe for each gene; best i.e. significant and largest abs(logFC)


pred_all_ann=left_join(pred_all,results_mRNA_short_sig,by=c("target_symbol"="SYMBOL"))

pred.names=pred_all_ann$target_symbol
pred.probes=pred_all_ann$PROBEID
pred.probes=unique(pred.probes) # just in case

exp_norm = Biobase::exprs(data_rma_norm)
exp_norm_hm=exp_norm[rownames(exp_norm)%in%pred.probes,]
exp_norm_hm= exp_norm_hm[match(pred.probes, row.names(exp_norm_hm)),,drop=FALSE]
exp_norm_hm_scaled=exp_norm_hm-rowMeans(exp_norm_hm)

#change probeIDS to gene IDS
sel=pred_all_ann[,c(3,5,7)][pred_all_ann$PROBEID%in%pred.probes,]

#this above preserves the order of elements
#make sure that the order of sel and probeids from rownames is the same
sel_reord=sel[match(pred.probes,sel$PROBEID),]
rownames(exp_norm_hm_scaled)=sel_reord$target_symbol


hm_mRNA=ComplexHeatmap::pheatmap(exp_norm_hm_scaled ,color = colour_palette_hm,
  fontsize=6, fontsize_row=7,cellheight=7.5,
  annotation_col=annot_col,annotation_colors=ann_colors)





## ---- table_toptargets

top_targets.mir.i=pred.cons.i_ann[,c(3,6,8,10,13)]
top_targets.mir.i=arrange(top_targets.mir.i, desc(abs(logFC)) )

table_toptargets=as.data.frame(top_targets.mir.i[1:5,])

caption.tab=paste("Top up-regulated targets (by predicted conserved sites) of microRNA",i,
  "(GSE12049 data set).", sep=" ")

rws=seq(1, (nrow(table_toptargets)-1), by = 2)
col=rep("\\rowcolor[gray]{0.95}", length(rws))

print(xtable(as.data.frame(table_toptargets) , align=c("p{0.1cm}","c","c","p{4cm}","c","c"), 
display=c("s","s","d","s","f","e"),  #digits=c(0,0,3,3,3),
caption=caption.tab),
include.rownames=FALSE,floating=TRUE,
table.placement="H", caption.placement = "bottom",
add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")








## ---- hm_DE_mRNA_mir1
print(hm_mRNA)


## ---- GO_terms_report

# pred.cons.i_ann

# > pred.cons.i_ann
# # A tibble: 36 x 13
# # Groups:   target_symbol, mature_mirna_id [31]

# > head(as.data.frame(pred.cons.i_ann),n=2)
#   mature_mirna_acc mature_mirna_id target_symbol target_entrez     target_ensembl predicted.sum
# 1     MIMAT0000123   mmu-miR-1a-3p         Abca1         11303 ENSMUSG00000015243             2
# 2     MIMAT0000123   mmu-miR-1a-3p         Actn2         11472 ENSMUSG00000052374             1
#      PROBEID                                            GENENAME            ENSEMBL     logFC
# 1 1421840_at ATP-binding cassette, sub-family A (ABC1), member 1 ENSMUSG00000015243 1.1972806
# 2 1448327_at                                     actinin alpha 2 ENSMUSG00000052374 0.9472434
#     AveExpr     P.Value  adj.P.Val
# 1  9.338973 0.002046557 0.08594378
# 2 11.575566 0.001449803 0.07665576

# > unique(pred.cons.i_ann$mature_mirna_id)
# [1] "mmu-miR-1a-3p"

#test set: DE targets of selected miRNA
DE_targets_ann_miR.i.probeIDs=unique(pred.cons.i_ann$PROBEID)


#match the background set
back_genes_idx <- genefilter::genefinder(data_final, 
                                        as.character(DE_targets_ann_miR.i.probeIDs), 
                                        method = "manhattan", scale = "none")

back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)

back_genes <- featureNames(data_final)[back_genes_idx]
# > length(back_genes)
# [1] 900
back_genes <- setdiff(back_genes, DE_targets_ann_miR.i.probeIDs)

# #plot with mean expression on the x-axis and curves for all genes, foreground genes and background genes, respectively
# library(geneplotter)
# multidensity(list(
#         all = results_mRNA[,"AveExpr"] ,
#         fore = results_mRNA[DE_targets_ann_miR.i.probeIDs , "AveExpr"],
#         back = results_mRNA[rownames(results_mRNA) %in% back_genes, "AveExpr"]),
#         col = c("#e46981", "#ae7ee2", "#a7ad4a"),
#      xlab = "mean expression",
#    main = "DE genes for CD-background-matching")


all_gene_IDs <- rownames(results_mRNA)
in_universe <- all_gene_IDs %in% c(DE_targets_ann_miR.i.probeIDs, back_genes)
in_selection <- all_gene_IDs %in% DE_targets_ann_miR.i.probeIDs

all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- all_gene_IDs[in_universe] 

top_GO_data = new("topGOdata", ontology = "BP", allGenes = all_genes,
 nodeSize = 10, annot = annFUN.db, affyLib = "mouse4302.db")


# For the method classic each GO category is tested independently
result_top_GO_classic = runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

nSigGO=sum(score(result_top_GO_classic)<=0.01, na.rm=TRUE)

#categories, top 54 (number of sig categories)
res_top_GO = GenTable(top_GO_data, 
        Fisher.classic = result_top_GO_classic,
        orderBy = "Fisher.classic" , topNodes = nSigGO)


## ---- table_topGO

table_topGO=as.data.frame(res_top_GO[1:5,c(1,2,6)])

caption.tab=paste("Top GO terms overrepresented in up-regulated targets (by predicted conserved sites) of microRNA",i,
  "(GSE12049 data set). P value is from one sided Fisher test.", sep=" ")

rws=seq(1, (nrow(table_topGO)-1), by = 2)
col=rep("\\rowcolor[gray]{0.95}", length(rws))

print(xtable(as.data.frame(table_topGO) , align=c("p{0.1cm}","c","p{7cm}","c"), 
display=c("s","s","s","e"),  #digits=c(0,0,3,3,3),
caption=caption.tab),
include.rownames=FALSE,floating=TRUE,
table.placement="H", caption.placement = "bottom",
add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")

## ---- plot_GO
showSigOfNodes(top_GO_data, score(result_top_GO_classic), firstSigNodes = 3, useInfo = 'def')

#to save to file
#fnamepth=file.path(resdir,"DAG.GO.pdf")
#printGraph(top_GO_data, result_top_GO_classic, firstSigNodes = 5, fn.prefix = fnamepth, useInfo = "all", pdfSW = TRUE)

## ---- reactome_enrich
entrez_ids <- mapIds(mouse4302.db, 
      keys = rownames(results_mRNA), 
      keytype = "PROBEID",
      column = "ENTREZID")



reactome_enrich <- enrichPathway(gene = entrez_ids[DE_targets_ann_miR.i.probeIDs], 
                                universe = entrez_ids[c(DE_targets_ann_miR.i.probeIDs, 
                                                        back_genes)],
                                organism = "mouse",
                                pvalueCutoff = 0.1,
                                qvalueCutoff = 0.9, 
                                readable = TRUE)

res_top_react=as.data.frame(reactome_enrich)[1:6]
res_top_react=res_top_react[res_top_react$p.adjust<0.1,]


## ---- table_reactome
react_table_rep=res_top_react[1:5,c(1,2,5,6)]

caption.tab=paste("Overrpresented Reactome pathways in up-regulated targets of ",i,". P value is from one sided Fisher test.", sep="")

rws=seq(1, (nrow(react_table_rep)-1), by = 2)
col=rep("\\rowcolor[gray]{0.95}", length(rws))

print(xtable(as.data.frame(react_table_rep) , align=c("p{0.1cm}","c","p{7cm}","c","c"), 
display=c("s","s","s","e","e"),  digits=c(0,0,0,3,3),
caption=caption.tab),
include.rownames=FALSE,floating=TRUE,
table.placement="H", caption.placement = "bottom",
add.to.row = list(pos = as.list(rws), command = col), size="\\scriptsize")


## ---- plot_reactome
reactome_enrich <- pairwise_termsim(reactome_enrich)
emapplot(reactome_enrich, showCategory = 20)





## the following chunks produce the above output for all significant miRNAs


## ----- sel_miRNA_dirs

dir_miRNAsel_pth=file.path(resdir,file.path(mi_mRNA_int.dirname,"significant-miRNAs"))
dir.create(dir_miRNAsel_pth)


mirDN_mrnaUP_dir=file.path(dir_miRNAsel_pth,"miR_dn.mRNA_up")

mirUP_mrnaDN_dir=file.path(dir_miRNAsel_pth,"miR_up.mRNA_dn")

mirDN_mrnaDN_dir=file.path(dir_miRNAsel_pth,"miR_dn.mRNA_dn")

mirUP_mrnaUP_dir=file.path(dir_miRNAsel_pth,"miR_up.mRNA_up")


#other objects are already prepared earlier for the examples in the report

# results_mRNA_short=results_mRNA[,c(1,2,3,4,5,6,8,9)]


# results_mRNA_short_sig=results_mRNA_short%>%
#   filter(adj.P.Val<FDR_CO)%>%
#   filter(abs(logFC)>lfc_CO)%>%
#   group_by(SYMBOL)%>%
#   slice(which.max(abs(logFC)))
# colour_palette_hm=colorRamp2(c(-2, 0, 2), c("dodgerblue4", "ivory", "chocolate3")) #z-scored

# annot_col=samples
# annot_col=as.data.frame(annot_col[,-1])
# colnames(annot_col)="Strain"

# ann_colors=list(
#   Strain=c(Mut="darkorange2", WT="dodgerblue4"))

## ---- selected_miRNA_analysis_functions




## ---- mirDN_mrnaUP_sig
dir_selected_group=mirDN_mrnaUP_dir
dir.create(dir_selected_group)

mir_to_process=unique(c(sig.miRdn_GENup.predCons$miRNA,
    sig.miRdn_GENup.predAll$miRNA,
    sig.miRdn_GENup.predCons.comb$miRNA,
    sig.miRdn_GENup.predAll.comb$miRNA))

pred.cons=miRdn_GENup_pred_cons%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

pred.all=miRdn_GENup_pred_all%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

valid=miRdn_GENup_valid

# select targets of each mir_to_process from here
GENup_pred_cons =miRdn_GENup_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_cons =miRdn_GENdn_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))

GENup_pred_all =miRdn_GENup_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_all =miRdn_GENdn_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))



for (i in mir_to_process){

  outdir.i=file.path(dir_selected_group,i)
  dir.create(outdir.i)

  #prep data
  pred.cons.i=pred.cons[pred.cons$mature_mirna_id==i,]
  pred.all.i=pred.all[pred.all$mature_mirna_id==i,]
  valid.i=valid[valid$mature_mirna_id==i,]

  pred.cons.i_ann=left_join(pred.cons.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  pred.all.i_ann=left_join(pred.all.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  valid.i_ann=left_join(valid.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  fname1=paste("predicted_conserved_targets",i,"tsv",sep=".")
  write.table(pred.cons.i_ann, file.path(outdir.i,fname1), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname2=paste("predicted_all_targets",i,"tsv",sep=".")
  write.table(pred.all.i_ann, file.path(outdir.i,fname2), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname3=paste("validated_targets",i,"tsv",sep=".")
  write.table(valid.i_ann, file.path(outdir.i,fname3), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


  #all targets lumped together, the best proble selected (lowest FDR)
  targets.i_all_ann=rbind(pred.cons.i_ann,pred.all.i_ann,valid.i_ann)
  targets.i_all_ann=targets.i_all_ann%>% group_by(target_symbol)%>% slice(which.min(adj.P.Val))

  expr_miRNA.i=as.numeric(exp_mir[rownames(exp_mir)==i,])


  #targets of miR.i amongst up and down regulated genes
  GENup_pred_cons.i=GENup_pred_cons[GENup_pred_cons$mature_mirna_id==i,]
  GENdn_pred_cons.i=GENdn_pred_cons[GENdn_pred_cons$mature_mirna_id==i,]
  GENup_pred_all.i=GENup_pred_all[GENup_pred_all$mature_mirna_id==i,]
  GENdn_pred_all.i=GENdn_pred_all[GENdn_pred_all$mature_mirna_id==i,]

  pred_cons=rbind(GENup_pred_cons.i,GENdn_pred_cons.i)
  pred_all=rbind(GENup_pred_all.i,GENdn_pred_all.i)


  # plots

  caption_levels_boxplot.x=paste("Median normalised log expression of\n",i,sep=" ")

  plot_levels_boxplot(exp_miRNA=expr_miRNA.i,exp_mRNA=data_final,targets=targets.i_all_ann,dirname=outdir.i,caption.x=caption_levels_boxplot.x)


  ## for the heatmap select best probe for each gene; best i.e. significant and largest abs(logFC)
  ## these are in results_mRNA_short_sig

 n_cons=nrow(pred_cons)
  if(n_cons>1){
    caption_hm.i=paste("Expression of conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"conserved_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_cons,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  n_all=nrow(pred_all)
  if(n_all>1){
    caption_hm.i=paste("Expression of conserved and non-conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"all_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_all,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  #GO terms
  de_probes=unique(pred.cons.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"conserved_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }

  de_probes=unique(pred.all.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"all_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }

}


## ---- mirUP_mrnaDN_sig

dir_selected_group=mirUP_mrnaDN_dir
dir.create(dir_selected_group)

mir_to_process=unique(c(sig.miRup_GENdn.predCons$miRNA,
    sig.miRup_GENdn.predAll$miRNA,
    sig.miRup_GENdn.predCons.comb$miRNA,
    sig.miRup_GENdn.predAll.comb$miRNA))

pred.cons=miRup_GENdn_pred_cons%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

pred.all=miRup_GENdn_pred_all%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

valid=miRup_GENdn_valid

# select targets of each mir_to_process from here
GENup_pred_cons =miRup_GENup_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_cons =miRup_GENdn_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))

GENup_pred_all =miRup_GENup_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_all =miRup_GENdn_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))



for (i in mir_to_process){

  outdir.i=file.path(dir_selected_group,i)
  dir.create(outdir.i)

  #prep data
  pred.cons.i=pred.cons[pred.cons$mature_mirna_id==i,]
  pred.all.i=pred.all[pred.all$mature_mirna_id==i,]
  valid.i=valid[valid$mature_mirna_id==i,]

  pred.cons.i_ann=left_join(pred.cons.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  pred.all.i_ann=left_join(pred.all.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  valid.i_ann=left_join(valid.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)


  fname1=paste("predicted_conserved_targets",i,"tsv",sep=".")
  write.table(pred.cons.i_ann, file.path(outdir.i,fname1), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname2=paste("predicted_all_targets",i,"tsv",sep=".")
  write.table(pred.all.i_ann, file.path(outdir.i,fname2), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname3=paste("validated_targets",i,"tsv",sep=".")
  write.table(valid.i_ann, file.path(outdir.i,fname3), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


  #all targets lumped together, the best proble selected (lowest FDR)
  targets.i_all_ann=rbind(pred.cons.i_ann,pred.all.i_ann,valid.i_ann)
  targets.i_all_ann=targets.i_all_ann%>% group_by(target_symbol)%>% slice(which.min(adj.P.Val))

  expr_miRNA.i=as.numeric(exp_mir[rownames(exp_mir)==i,])


  #targets of miR.i amongst up and down regulated genes
  GENup_pred_cons.i=GENup_pred_cons[GENup_pred_cons$mature_mirna_id==i,]
  GENdn_pred_cons.i=GENdn_pred_cons[GENdn_pred_cons$mature_mirna_id==i,]
  GENup_pred_all.i=GENup_pred_all[GENup_pred_all$mature_mirna_id==i,]
  GENdn_pred_all.i=GENdn_pred_all[GENdn_pred_all$mature_mirna_id==i,]

  pred_cons=rbind(GENup_pred_cons.i,GENdn_pred_cons.i)
  pred_all=rbind(GENup_pred_all.i,GENdn_pred_all.i)


  # plots

  caption_levels_boxplot.x=paste("Median normalised log expression of\n",i,sep=" ")

  plot_levels_boxplot(exp_miRNA=expr_miRNA.i,exp_mRNA=data_final,targets=targets.i_all_ann,dirname=outdir.i,caption.x=caption_levels_boxplot.x)


  ## for the heatmap select best probe for each gene; best i.e. significant and largest abs(logFC)
  ## these are in results_mRNA_short_sig

 n_cons=nrow(pred_cons)
  if(n_cons>1){
    caption_hm.i=paste("Expression of conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"conserved_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_cons,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  n_all=nrow(pred_all)
  if(n_all>1){
    caption_hm.i=paste("Expression of conserved and non-conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"all_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_all,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  #GO terms
  de_probes=unique(pred.cons.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"conserved_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }

  de_probes=unique(pred.all.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"all_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }
  

}


## ---- mirDN_mrnaDN_sig

dir_selected_group=mirDN_mrnaDN_dir
dir.create(dir_selected_group)

mir_to_process=unique(c(sig.miRdn_GENdn.predCons$miRNA,
    sig.miRdn_GENdn.predAll$miRNA,
    sig.miRdn_GENdn.predCons.comb$miRNA,
    sig.miRdn_GENdn.predAll.comb$miRNA))

pred.cons=miRdn_GENdn_pred_cons%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

pred.all=miRdn_GENdn_pred_all%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

valid=miRdn_GENdn_valid

# select targets of each mir_to_process from here
GENup_pred_cons =miRdn_GENup_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_cons =miRdn_GENdn_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))

GENup_pred_all =miRdn_GENup_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_all =miRdn_GENdn_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))



for (i in mir_to_process){

  outdir.i=file.path(dir_selected_group,i)
  dir.create(outdir.i)

  #prep data
  pred.cons.i=pred.cons[pred.cons$mature_mirna_id==i,]
  pred.all.i=pred.all[pred.all$mature_mirna_id==i,]
  valid.i=valid[valid$mature_mirna_id==i,]

  pred.cons.i_ann=left_join(pred.cons.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  pred.all.i_ann=left_join(pred.all.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  valid.i_ann=left_join(valid.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)


  fname1=paste("predicted_conserved_targets",i,"tsv",sep=".")
  write.table(pred.cons.i_ann, file.path(outdir.i,fname1), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname2=paste("predicted_all_targets",i,"tsv",sep=".")
  write.table(pred.all.i_ann, file.path(outdir.i,fname2), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname3=paste("validated_targets",i,"tsv",sep=".")
  write.table(valid.i_ann, file.path(outdir.i,fname3), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


  #all targets lumped together, the best proble selected (lowest FDR)
  targets.i_all_ann=rbind(pred.cons.i_ann,pred.all.i_ann,valid.i_ann)
  targets.i_all_ann=targets.i_all_ann%>% group_by(target_symbol)%>% slice(which.min(adj.P.Val))

  expr_miRNA.i=as.numeric(exp_mir[rownames(exp_mir)==i,])


  #targets of miR.i amongst up and down regulated genes
  GENup_pred_cons.i=GENup_pred_cons[GENup_pred_cons$mature_mirna_id==i,]
  GENdn_pred_cons.i=GENdn_pred_cons[GENdn_pred_cons$mature_mirna_id==i,]
  GENup_pred_all.i=GENup_pred_all[GENup_pred_all$mature_mirna_id==i,]
  GENdn_pred_all.i=GENdn_pred_all[GENdn_pred_all$mature_mirna_id==i,]

  pred_cons=rbind(GENup_pred_cons.i,GENdn_pred_cons.i)
  pred_all=rbind(GENup_pred_all.i,GENdn_pred_all.i)


  # plots

  caption_levels_boxplot.x=paste("Median normalised log expression of\n",i,sep=" ")

  plot_levels_boxplot(exp_miRNA=expr_miRNA.i,exp_mRNA=data_final,targets=targets.i_all_ann,dirname=outdir.i,caption.x=caption_levels_boxplot.x)


  ## for the heatmap select best probe for each gene; best i.e. significant and largest abs(logFC)
  ## these are in results_mRNA_short_sig

  n_cons=nrow(pred_cons)
  if(n_cons>1){
    caption_hm.i=paste("Expression of conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"conserved_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_cons,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  n_all=nrow(pred_all)
  if(n_all>1){
    caption_hm.i=paste("Expression of conserved and non-conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"all_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_all,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  #GO terms
  de_probes=unique(pred.cons.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"conserved_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }

  de_probes=unique(pred.all.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"all_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }

}

## ---- mirUP_mrnaUP_sig

dir_selected_group=mirUP_mrnaUP_dir
dir.create(dir_selected_group)

mir_to_process=unique(c(sig.miRup_GENup.predCons$miRNA,
    sig.miRup_GENup.predAll$miRNA,
    sig.miRup_GENup.predCons.comb$miRNA,
    sig.miRup_GENup.predAll.comb$miRNA))

pred.cons=miRup_GENup_pred_cons%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

pred.all=miRup_GENup_pred_all%>% 
  group_by(target_symbol,mature_mirna_id)%>% 
  slice(which.max(predicted.sum))

valid=miRup_GENup_valid

# select targets of each mir_to_process from here
GENup_pred_cons =miRup_GENup_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_cons =miRup_GENdn_pred_cons%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))

GENup_pred_all =miRup_GENup_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))
GENdn_pred_all =miRup_GENdn_pred_all%>% group_by(target_symbol,mature_mirna_id)%>% slice(which.max(predicted.sum))



for (i in mir_to_process){

  outdir.i=file.path(dir_selected_group,i)
  dir.create(outdir.i)

  #prep data
  pred.cons.i=pred.cons[pred.cons$mature_mirna_id==i,]
  pred.all.i=pred.all[pred.all$mature_mirna_id==i,]
  valid.i=valid[valid$mature_mirna_id==i,]

  pred.cons.i_ann=left_join(pred.cons.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  pred.all.i_ann=left_join(pred.all.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  valid.i_ann=left_join(valid.i,results_mRNA_short,by=c("target_symbol"="SYMBOL")) %>%
    filter(adj.P.Val<FDR_CO) %>%
    filter(abs(logFC)>lfc_CO)

  fname1=paste("predicted_conserved_targets",i,"tsv",sep=".")
  write.table(pred.cons.i_ann, file.path(outdir.i,fname1), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname2=paste("predicted_all_targets",i,"tsv",sep=".")
  write.table(pred.all.i_ann, file.path(outdir.i,fname2), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  fname3=paste("validated_targets",i,"tsv",sep=".")
  write.table(valid.i_ann, file.path(outdir.i,fname3), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)


  #all targets lumped together, the best proble selected (lowest FDR)
  targets.i_all_ann=rbind(pred.cons.i_ann,pred.all.i_ann,valid.i_ann)
  targets.i_all_ann=targets.i_all_ann%>% group_by(target_symbol)%>% slice(which.min(adj.P.Val))

  expr_miRNA.i=as.numeric(exp_mir[rownames(exp_mir)==i,])


  #targets of miR.i amongst up and down regulated genes
  GENup_pred_cons.i=GENup_pred_cons[GENup_pred_cons$mature_mirna_id==i,]
  GENdn_pred_cons.i=GENdn_pred_cons[GENdn_pred_cons$mature_mirna_id==i,]
  GENup_pred_all.i=GENup_pred_all[GENup_pred_all$mature_mirna_id==i,]
  GENdn_pred_all.i=GENdn_pred_all[GENdn_pred_all$mature_mirna_id==i,]

  pred_cons=rbind(GENup_pred_cons.i,GENdn_pred_cons.i)
  pred_all=rbind(GENup_pred_all.i,GENdn_pred_all.i)


  # plots

  caption_levels_boxplot.x=paste("Median normalised log expression of\n",i,sep=" ")

  plot_levels_boxplot(exp_miRNA=expr_miRNA.i,exp_mRNA=data_final,targets=targets.i_all_ann,dirname=outdir.i,caption.x=caption_levels_boxplot.x)


  ## for the heatmap select best probe for each gene; best i.e. significant and largest abs(logFC)
  ## these are in results_mRNA_short_sig

  n_cons=nrow(pred_cons)
  if(n_cons>1){
    caption_hm.i=paste("Expression of conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"conserved_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_cons,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  n_all=nrow(pred_all)
  if(n_all>1){
    caption_hm.i=paste("Expression of conserved and non-conserved DE targets of ",i,".\nNormalised log expression is plotted after scaling by row.", sep="")
    hm_fname=paste("heatmap_exprs",i,"all_targets","pdf", sep=".")

    plot_heatmap(genes_to_plot=pred_all,res_mRNA=results_mRNA_short_sig,expr_set=data_rma_norm,fname=hm_fname,dirname=outdir.i,caption=caption_hm.i)
  }

  #GO terms
  de_probes=unique(pred.cons.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"conserved_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }

  de_probes=unique(pred.all.i_ann$PROBEID)
  nprobes=length(de_probes)
  if(nprobes>1){
    fname_GO=paste(i,"all_targets",sep=".")

    get_GO_terms_react(de_probes=de_probes,expr_set=data_final,results_mRNA=results_mRNA,fname=fname_GO, dirname=outdir.i)
  }

}


