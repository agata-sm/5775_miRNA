plot_levels_boxplot <- function(exp_miRNA,exp_mRNA,targets,dirname,caption.x){

  miRNA_median.wt=median(exp_miRNA[6:10])
  miRNA_median.mut=median(exp_miRNA[1:5])
  
  targets.to_plot_symbol=targets$target_symbol
  targets.to_plot_probes=targets$PROBEID
  names(targets.to_plot_probes)=targets.to_plot_symbol

  outdir_boxplots=file.path(dirname,paste("boxplots-levels",i, sep="."))
  dir.create(outdir_boxplots)


  n_probs=length(targets.to_plot_probes)

  for (jj in (1:n_probs)){

    j=targets.to_plot_probes[jj]
    name.j=names(targets.to_plot_probes[jj])

    caption_levels.y=paste("Normalised log expression of\n",name.j,sep=" ")

    mRNA_exprs.j=exprs(data_final[rownames(data_final)==j,])

    df_plot.j=as.data.frame(t(mRNA_exprs.j))
    colnames(df_plot.j)="mRNA"

    df_plot.j$strain=c(rep("WT",3),rep("Mut",3))

    df_plot.j$miRNA_median=c(rep(miRNA_median.wt,3),rep(miRNA_median.mut,3))

    caption.j=paste("Expression of",name.j,"vs. median expression of",i,"in WT and Mut.")

    plot.j_save=ggplot(df_plot.j, aes(x=miRNA_median, y=mRNA,colour=strain)) + geom_boxplot(width=0.5) + geom_point() +
    scale_color_manual(values = c("darkorange2", "dodgerblue4")) + theme_pca +
    xlab(caption.x)+ ylab(caption_levels.y) + ggtitle(caption.j)

    fpath=(file.path(outdir_boxplots, paste(name.j,i,"pdf",sep=".")))
    pdf(fpath)
    print(plot.j_save)
    dev.off()

  }

}



plot_heatmap <-function(genes_to_plot,res_mRNA,expr_set,dirname,fname,caption){

  #get probe IDS and names
  genes_all_ann=left_join(genes_to_plot,res_mRNA,by=c("target_symbol"="SYMBOL"))

  genes.names=genes_all_ann$target_symbol
  genes.probes=genes_all_ann$PROBEID
  genes.probes=unique(genes.probes) # just in case

  #get normalised expression values
  exp_norm = Biobase::exprs(expr_set)
  exp_norm_hm=exp_norm[rownames(exp_norm)%in%genes.probes,]
  exp_norm_hm= exp_norm_hm[match(genes.probes, row.names(exp_norm_hm)),,drop=FALSE]
  exp_norm_hm_scaled=exp_norm_hm-rowMeans(exp_norm_hm)

  #change probeIDS to gene IDS
  sel=genes_all_ann[,c(3,5,7)][genes_all_ann$PROBEID%in%genes.probes,]
  #this above preserves the order of elements
  #make sure that the order of sel and probeids from rownames is the same
  sel_reord=sel[match(genes.probes,sel$PROBEID),]
  rownames(exp_norm_hm_scaled)=sel_reord$target_symbol

  hm_mRNA=ComplexHeatmap::pheatmap(exp_norm_hm_scaled ,color = colour_palette_hm,
  fontsize=6, fontsize_row=7,cellheight=7.5,
  annotation_col=annot_col,annotation_colors=ann_colors,
  main=caption)

  fpath=file.path(dirname,fname)
  pdf(fpath)
  print(hm_mRNA)
  dev.off()

}


get_GO_terms_react <-function(de_probes, expr_set,results_mRNA, fname, dirname){

  #match the background set
  back_genes_idx <- genefilter::genefinder(expr_set, as.character(de_probes), 
                                        method = "manhattan", scale = "none")
  back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
  back_genes <- featureNames(expr_set)[back_genes_idx]
  back_genes <- setdiff(back_genes, de_probes)

  all_gene_IDs <- rownames(results_mRNA)
  in_universe <- all_gene_IDs %in% c(de_probes, back_genes)
  in_selection <- all_gene_IDs %in% de_probes

  all_genes <- in_selection[in_universe]
  all_genes <- factor(as.integer(in_selection[in_universe]))
  names(all_genes) <- all_gene_IDs[in_universe] 

  top_GO_data = new("topGOdata", ontology = "BP", allGenes = all_genes,
  nodeSize = 10, annot = annFUN.db, affyLib = "mouse4302.db")


  result_top_GO_classic = runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")

  nSigGO=sum(score(result_top_GO_classic)<=0.01, na.rm=TRUE)

  res_top_GO = GenTable(top_GO_data, 
        Fisher.classic = result_top_GO_classic,
        orderBy = "Fisher.classic" , topNodes = nSigGO)

  fname_GOtab=paste("GO-terms",fname,"tsv",sep=".")
  write.table(res_top_GO, file.path(dirname,fname_GOtab), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  if (nSigGO >0){
    fname_plot_GO=file.path(dirname,paste("GO-terms",fname,"DAG","pdf",sep="."))

    printGraph(top_GO_data, result_top_GO_classic, firstSigNodes = 5, fn.prefix = fname_plot_GO, useInfo = "all", pdfSW = TRUE)
  }
  entrez_ids <- mapIds(mouse4302.db, 
      keys = rownames(results_mRNA), 
      keytype = "PROBEID",
      column = "ENTREZID")

  
  reactome_enrich <- enrichPathway(gene = entrez_ids[de_probes], 
                                universe = entrez_ids[c(de_probes, 
                                                        back_genes)],
                                organism = "mouse",
                                pvalueCutoff = 0.1,
                                qvalueCutoff = 0.9, 
                                readable = TRUE)

  res_top_react=as.data.frame(reactome_enrich)
  res_top_react=res_top_react[res_top_react$p.adjust<0.1,]

  fname_Reacttab=paste("Reactome",fname,"tsv",sep=".")
  write.table(res_top_react, file.path(dirname,fname_Reacttab), sep="\t",col.names=TRUE, row.names=FALSE, quote=FALSE)

  rws=nrow(res_top_react)
  if (rws>1){
    reactome_enrich <- pairwise_termsim(reactome_enrich)

    fname_plot_react=file.path(dirname,paste("Reactome",fname,"emap","pdf",sep="."))
    pdf(fname_plot_react)
    print(emapplot(reactome_enrich, showCategory = 20))
    dev.off()
  }
}






