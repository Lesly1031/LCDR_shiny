limma_choices <- function(phenotype){
  phenotype = data.frame(phenotype)
  phenotype1 = phenotype[, sapply(phenotype, function(col) length(unique(col)) > 1), drop = FALSE]
  
  if(ncol(phenotype1) == 0) {
    warning("No columns left in phenotype1 after filtering.")
    return(NULL)  # Exit if no columns meet the criterion
  }
  
  choices_comb = lapply(seq_len(ncol(phenotype1)), function(x){
    # Now use the column name directly for grouping
    res = dplyr::count(phenotype1, !!sym(colnames(phenotype1)[x]))
    if(ncol(res) == 2) {
      names(res) = c("name", "n")
    } else {
      warning("Expected 'res' to have 2 columns, found: ", ncol(res))
      return(NA)
    }
    choices = res$name[which(res$n >= 2)]
    
    tryCatch({
      combin = matrix(combn(choices, 2), nrow = 2)
      contrast1 = sapply(1:ncol(combin), function(i){
        paste0(combin[,i][1], " vs ", combin[,i][2])
      })
      contrast2 = sapply(1:ncol(combin), function(i){
        paste0(combin[,i][2], " vs ", combin[,i][1])
      })
      
      c(paste0(colnames(phenotype1)[x], ":", contrast1), paste0(colnames(phenotype1)[x], ":", contrast2))
    }, error = function(e){
      warning("Error processing combinations: ", e$message)
      NA
    })
  })
  
  choices_comb = choices_comb[!is.na(choices_comb)]
  final_choices = unlist(choices_comb)
  return(final_choices)
}



limma_table = function(input,output,session, counts,pheno,group_info){
  # group_info = "Group:group1 vs group2"
  # counts = read_fst("M:/dept/Dept_BBSR/Projects/Fridley_Brooke/1959_LCDR/work/chip_seq_shiny/chip_seq/files/study/demo/bw_read_counts.fst")
  # pheno = fread("M:/dept/Dept_BBSR/Projects/Fridley_Brooke/1959_LCDR/work/chip_seq_shiny/chip_seq/files/study/demo/phenotype.csv")
 # counts = read_fst("M:/dept/Dept_BBSR/Projects/Fridley_Brooke/1959_LCDR/work/chip_seq_shiny/chip_seq/files/study/3675/bw_read_counts.fst")
#  pheno = fread("M:/dept/Dept_BBSR/Projects/Fridley_Brooke/1959_LCDR/work/chip_seq_shiny/chip_seq/files/study/3675/phenotype.csv")
  
  counts = data.frame(counts,check.names=F)
  pheno = data.frame(pheno,check.names=F)
 # counts = data.frame(na.omit(counts))
  
  pheno$X[pheno$X==""]=NA
  pheno = pheno[complete.cases(pheno$X),]
  counts = data.frame(counts[,c(1:6,which(names(counts)%in% pheno$X))],check.names = F)
  
  counts_res = counts[,-c(1:6)]
  counts_res_div = lapply(1:length(counts_res),function(x){
    counts <- round(counts_res[[x]] / 70, 0) # ??????? am I correct????
  })
  
  ## Explore a little portion of the count matrix 
  counts_res_div = do.call(cbind,counts_res_div)
  rownames(counts_res_div) = counts$GeneID

  
  dge <- DGEList(counts=counts_res_div)
  
  #group = droplevels(pheno$Group)
  group_info1 = str_split(input[[group_info]],":")[[1]][1]
  comp_group = str_split(input[[group_info]],":")[[1]][2]
  comp_group_sub = make.names(str_split(comp_group," vs ")[[1]])
  
  grp  <- factor(make.names(pheno[[group_info1]]))
  
  design <- model.matrix(~ 0+grp)
  colnames(design) = make.names(colnames(design))
  
  keep <- filterByExpr(dge, design)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  v <- voom(dge,design,plot=F)
  contr <- makeContrasts(contrasts = paste0(paste0("grp",comp_group_sub[1]),"-",paste0("grp",comp_group_sub[2])), levels = design)
  #fit1 <- lmFit(dge$counts,design)
  fit1 <- lmFit(v,design,weights=v$weights)
  fit2 <- contrasts.fit(fit1, contr)
  fitW <- eBayes(fit2)
  res=data.frame(fitW)
  res=data.frame(res[,c("coefficients","Amean","t","lods","p.value")])
  res$padj=p.adjust(res$p.value,method="BH")
  names(res)=c("log2FoldChange","baseMean","t stat","log odds of differential expression","pvalue","padj")
  
  rownames(res) = rownames(dge$counts)
  res = merge(counts[,1:5],res,by.y="row.names",by.x="GeneID")
  for(i in 7:ncol(res)){
    res[,i] = round(res[,i],4)
  }
  res=res[order(res$pvalue ,res$`log odds of differential expression`,res$GeneName,decreasing = F),]
  
  res
}