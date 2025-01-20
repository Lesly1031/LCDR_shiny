library(tidyverse)
#  contrasts_func<-function(input,output,session,phenotype,group_info,comp_group){
# observeEvent(input[[group_info]],{
#   tryCatch({
#     phenotype=phenotype
#     choices <-  factor(phenotype[,colnames(phenotype)[colnames(phenotype)==input[[paste0("group_limma",x)]]]])
#     combin <- matrix(combn(levels(choices), 2),nrow=2)
#     #combin <- (matrix(combn(levels(limma_var_pre()[[x]]), 2),nrow=2))
#     # combin=sapply(combin, function(xx) as.character(xx), simplify = FALSE)
#     contrast1=sapply(1:ncol(combin),function(i){
#       paste0(combin[,i][1],"-",combin[,i][2])
#     })
#     
#     contrast2=sapply(1:ncol(combin),function(i){
#       paste0(combin[,i][2],"-",combin[,i][1])
#     })
#     contrasts = c(contrast1,contrast2)
#     updateSelectizeInput(session, input[[comp_group]], choices = contrasts,server=T)
#   }, error=function(e){
#     "error"
#   })
# })
# 
#  }


library(tidyverse)
#  contrasts_func<-function(input,output,session,phenotype,group_info,comp_group){
# observeEvent(input[[group_info]],{
#   tryCatch({
#     phenotype=phenotype
#     choices <-  factor(phenotype[,colnames(phenotype)[colnames(phenotype)==input[[paste0("group_limma",x)]]]])
#     combin <- matrix(combn(levels(choices), 2),nrow=2)
#     #combin <- (matrix(combn(levels(limma_var_pre()[[x]]), 2),nrow=2))
#     # combin=sapply(combin, function(xx) as.character(xx), simplify = FALSE)
#     contrast1=sapply(1:ncol(combin),function(i){
#       paste0(combin[,i][1],"-",combin[,i][2])
#     })
#     
#     contrast2=sapply(1:ncol(combin),function(i){
#       paste0(combin[,i][2],"-",combin[,i][1])
#     })
#     contrasts = c(contrast1,contrast2)
#     updateSelectizeInput(session, input[[comp_group]], choices = contrasts,server=T)
#   }, error=function(e){
#     "error"
#   })
# })
# 
#  }

limma_choices <- function(phenotype){
  # exclude continuous variables
  phenotype = data.frame(phenotype)
  which_conti = sapply(names(phenotype),function(x) !is.numeric(phenotype[[x]]) )
  phenotype = data.frame(phenotype[names(which(which_conti))])
  
  phenotype1 =  phenotype[, sapply(phenotype, function(col) length(unique(col))) > 1]
  
  
  
  choices_comb = lapply(2:ncol(phenotype1),function(x){
    #choices=levels(factor(phenotype1[[colnames(phenotype1)[x]]]))
    res =  dplyr::count(phenotype1,phenotype1[[colnames(phenotype1)[x]]])
    res[res==""] = NA
    res = na.omit(res)
    names(res) = c("name","n")
    choices = res$name[which(res$n>=3)]
    
    tryCatch(      {combin <- (matrix(combn(choices, 2),nrow=2))
    contrast1=sapply(1:ncol(combin),function(i){
      paste0(combin[,i][1]," vs ",combin[,i][2])
    })
    
    contrast2=sapply(1:ncol(combin),function(i){
      paste0(combin[,i][2]," vs ",combin[,i][1])
    })
    c(paste0(colnames(phenotype1)[x],":",contrast1),paste0(colnames(phenotype1)[x],":",contrast2))
    },error = function(e){
      NA
    })
    
    
  })
  choices_comb = choices_comb[!is.na(choices_comb)]
  unlist(choices_comb)
}


limma_UI<-function(headings,output_info,download_info,group_info,phenotype,study_type){
    column(10,
           panel(
             status="info",heading=paste0("Study: ",headings),
             dropdown(
               selectizeInput(group_info, label="Select group annotation", 
                              choices =limma_choices(phenotype)),
               downloadButton(download_info, "Download"),
               circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
               tooltip = tooltipOptions(title = "Click to see inputs !")
             ),
             conditionalPanel( 
               condition="input.update_DEG==0",
               br(),
               p("Please click 'Run' button to generate results")
             ),
             
             conditionalPanel( 
               condition="input.update_DEG>0",
               br(),
               withLoader( DTOutput(output_info),type="html",loader="loader5")
               
             )
             #DTOutput(output_info)
             
           )
    )    
}



limma_server<-function(input,output,session,data,data_class,group_info,study_type){
  
  # observe({
  #   tryCatch({
  #     phenotype=data.frame(data_class)
  #     choices <-  factor(phenotype[,colnames(phenotype)[colnames(phenotype)==group_info]])
  #     combin <- matrix(combn(levels(choices), 2),nrow=2)
  #     #combin <- (matrix(combn(levels(limma_var_pre()[[x]]), 2),nrow=2))
  #     # combin=sapply(combin, function(xx) as.character(xx), simplify = FALSE)
  #     contrast1=sapply(1:ncol(combin),function(i){
  #       paste0(combin[,i][1],"-",combin[,i][2])
  #     })
  #
  #     contrast2=sapply(1:ncol(combin),function(i){
  #       paste0(combin[,i][2],"-",combin[,i][1])
  #     })
  #     contrasts = c(contrast1,contrast2)
  #     updateSelectizeInput(session, input[[comp_group]], choices = contrasts,server=T)
  #   }, error=function(e){
  #     "error"
  #   })
  # })
  # if(study_type=="metabolomics"){
  #   dt=data.frame(data)
  #   ID_new <- sapply( str_split(dt$ID, "_"), "[", 3)
  #   dt = dt[which(!is.na(ID_new)),]
  # }else if(grepl("RNAseq", study_type, fixed = TRUE)){
  #   dt=data.frame(data)
  #   ID_new <- sapply( str_split(dt$ID, "_"), "[", 2)
  #   dt = dt[which(!is.na(ID_new)),]
  #
  # }
  
  
  #tryCatch({

  dt=data.frame(data,check.names = F)
  

  if(grepl("Metabolomics",study_type)|grepl("Proteomics",study_type)|grepl("Lipidomics",study_type)){ ########### metabolomics ##################
   
    #ID_new <- sapply( str_split(dt$ID, "_"), "[", 1)
    ID_new = sapply(1:length(dt$ID),function(x){
      str_split(dt$ID[x], "_")[[1]][1]
    })
    dt = dt[which(!is.na(ID_new)),]
    
    dt_class=data.frame(data_class)
    
    dt_class$X[dt_class$X==""]=NA
    dt_class = dt_class[complete.cases(dt_class$X),]
    dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names = F)
    
    rownames_dt_class=as.character(dt_class[,1])
    colnames_dt_class=colnames(dt_class)[-1]
    dt_class=data.frame(dt_class[,-1])
    rownames(dt_class)=rownames_dt_class
    names(dt_class)=colnames_dt_class
    
    #group_info1 = str_split("experimental_group:Vehicle vs SR13800",":")[[1]][1]
    #comp_group = str_split("experimental_group:Vehicle vs SR13800",":")[[1]][2]
    
    group_info1 = str_split(input[[group_info]],":")[[1]][1]
    comp_group = str_split(input[[group_info]],":")[[1]][2]
    comp_group_sub = str_split(comp_group," vs ")[[1]]
    names_na=rownames(dt_class)[which(is.na(dt_class[,group_info1]))]
    dt_sub=dt
    rownames(dt_sub)=dt_sub$ID
    dt_sub=data.frame(dt_sub[,-1],check.names = F)
    
    if(length(names_na)!=0){
      dt_sub=dt_sub[,-(which(colnames(dt_sub)%in%names_na))]
      dt_class=na.omit(dt_class)
    } else{
      dt_sub=dt_sub
      dt_class=dt_class
    }
    #names(dt_sub)=make.names(names(dt_sub))
    dt_sub = data.frame(t(dt_sub),check.names=F) #### problematic ???
    #rownames(dt_class)= make.names(rownames(dt_class))
    dt_cb = merge(dt_sub,dt_class,by="row.names")
    dt_cb = subset(dt_cb,dt_cb[,group_info1]%in%comp_group_sub)
    dt_cb[[make.names(group_info1)]] = factor(dt_cb[[make.names(group_info1)]],levels =comp_group_sub )
    levels_class = levels(dt_cb[[make.names(group_info1)]])
    
    res1 =lapply(2:(ncol(dt_sub)+1),function(x){
      tryCatch({
        res=t.test(dt_cb[,x]~dt_cb[,make.names(group_info1)],var.equal =FALSE)
        c(colnames(dt_cb)[x],as.numeric(res$estimate),res$stderr,res$statistic,res$p.value)
      },error = function(e){
        NULL
      })
      
    })
    
    res = data.frame(do.call(rbind,res1),check.names = F)
    
    rownames(res) = res$V1;res = data.frame(res[,-1],check.names = F)
    
    for(i in 1:ncol(res)){
      res[,i] = as.numeric(res[,i])
    }
    names(res)= c(paste0("Estimated mean of ",levels_class[[1]]),
                  paste0("Estimated mean of ",levels_class[[2]]),"stderr","t_statistics","pval")
    
    res$padj =p.adjust(as.numeric(res$pval),"fdr")
    #res$baseMean = ((res[,paste0("Estimated mean of ",levels_class[[1]])])-(res[,paste0("Estimated mean of ",levels_class[[2]])]))
    
    fc = (res[,paste0("Estimated mean of ",levels_class[[1]])]) - (res[,paste0("Estimated mean of ",levels_class[[2]])])
    # fc = (res[,paste0("Estimated mean of ",comp_group_sub[[1]])]) - (res[,paste0("Estimated mean of ",comp_group_sub[[2]])])
    res$log2FoldChange = fc
    #res$baseMean = (1/2)*(log2(res[,paste0("Estimated mean of ",comp_group_sub[[1]])])+log2(res[,paste0("Estimated mean of ",comp_group_sub[[2]])]))
    res$baseMean = (1/2)*((res[,paste0("Estimated mean of ",levels_class[[1]])])+(res[,paste0("Estimated mean of ",levels_class[[2]])]))
    
    # res = res[,c(paste0("Estimated mean of ",comp_group_sub[[1]]),
    #              paste0("Estimated mean of ",comp_group_sub[[2]]),
    #              "baseMean","log2FoldChange",
    #              "stderr","t_statistics","pval","padj")
    # ]
    res = res[,c("baseMean","log2FoldChange",
                 "stderr","t_statistics","pval","padj")]
    names(res) = c("baseMean","log2FoldChange","stderr","t stat","pvalue","padj")
    res=res[order(res$pvalue,res$log2FoldChange,decreasing = F),]
    for(i in 1:4){
      res[,i] = data.frame(round(res[,i],4))
    }
    # split rownames
    res$rowID = rownames(res)
    res
  }else{ ########### RNAseq  ################
    ######### data process ################
  
    #ID_new <- sapply( str_split(dt$ID, "_"), "[", 1)
    ID_new = dt$ID
    dt = dt[which(!is.na(ID_new)),]
    
    dt_class=data.frame(data_class)
    
    dt_class$X[dt_class$X==""]=NA
    dt_class = dt_class[complete.cases(dt_class$X),]
    dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names = F)
    
    dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names=F)
    
    rownames_dt_class=as.character(dt_class[,1])
    colnames_dt_class=colnames(dt_class)[-1]
    dt_class=data.frame(dt_class[,-1])
    rownames(dt_class)=rownames_dt_class
    names(dt_class)=colnames_dt_class
    
    
    group_info1 = str_split(input[[group_info]],":")[[1]][1]
    comp_group = str_split(input[[group_info]],":")[[1]][2]
    comp_group_sub = str_split(comp_group," vs ")[[1]]
    # comp_group_sub = make.names(comp_group_sub)
    names_na=rownames(dt_class)[which(is.na(dt_class[,group_info1]))]
    dt_sub=dt
    rownames(dt_sub)=dt_sub$ID
    dt_sub=data.frame(dt_sub[,-1],check.names=F)
    
    
    
    
    if(length(names_na)!=0){
      dt_sub=dt_sub[,-(which(colnames(dt_sub)%in%names_na))]
      dt_class=na.omit(dt_class)
    } else{
      dt_sub=dt_sub
      dt_class=dt_class
    }
    #names(dt_sub)=make.names(names(dt_sub))
    
    
    #dt_class1=dt_class[colnames(dt_sub),]
    
    # for(i in 1:ncol(dt_sub)){
    #   dt_sub[,i]=log2(dt_sub[,i]+1)
    # }
    
    dt_class = dt_class[match(names(dt_sub),rownames(dt_class)),]
    #dt_class = data.frame(dt_class[match(names(dt_sub),rownames(dt_class)),])
    #names(dt_class) = group_info1
    
    x = tibble::rownames_to_column(dt_sub, "ID")
    
    #x <- DGEList(dt_sub) 
    grp  <- factor(dt_class[,group_info1])
    design <- model.matrix(~ 0+grp)
    colnames(design) = make.names(colnames(design))
    
    
    #keep.exprs <- filterByExpr(x, group = design)
    # keep.exprs<-apply(dt_sub,1,sum)>=10 & rowSums(dt_sub!=0)>=4
    # 
    # x <- x[keep.exprs, , keep.lib.sizes = FALSE]
    # x <- calcNormFactors(x,method="TMM") # calculate normalization factor: according to Chia-Ho's email on 01/11/2023
    # rnaseq_voomed = voom(x, design) 
    contr <- makeContrasts(contrasts = paste0(make.names(paste0("grp",comp_group_sub[1])),"-",make.names(paste0("grp",comp_group_sub[2]))), levels = design)
    fit1 = lmFit(dt_sub, design)
    fit2 <- contrasts.fit(fit1, contr)
    fitW <- eBayes(fit2)
    res=data.frame(fitW,check.names=F)
    res=data.frame(res[,c("coefficients","Amean","t","lods","p.value")],check.names=F)
    res$padj=p.adjust(res$p.value,method="BH")
    
    # colnames(design) = make.names(colnames(design))
    # y <- voom(na.omit(dt_sub),design,plot=F)
    # fit <- lmFit(y,design)
    # #  fit = lmFit(dt_sub,design)
    # 
    # #colnames(design) <- c(paste0("c",input$experiment_group),paste0("c",input$control_group))
    # #contr <- makeContrasts(contrasts = paste0(paste0("grp",str_split(comp_group[[1]],"-")[[1]][1]),"-",paste0("grp",str_split(comp_group[[1]],"-")[[1]][2])), levels = design)
    # 
    # ## With & without gene filtration
    # fit2 <- contrasts.fit(fit, contr)
    # fitW <- eBayes(fit2)
    # res=data.frame(fitW)
    # res=data.frame(res[,c("coefficients","Amean","t","lods","p.value")])
    # res$padj=p.adjust(res$p.value,method="BH")
    rownames(res) = rownames(fitW$coefficients)
    
    # res <- topTable(fitW,number=nrow(dt))
    #res$rowID = rownames(res)
    names(res)=c("log2FoldChange","baseMean","t stat","log odds of differential expression","pvalue","padj")
    res=res[order(res$pvalue,res$`log odds of differential expression`,decreasing = F),]
    #res = data.frame(res[,c("rowID","log2FoldChange","baseMean","t stat","log odds of differential expression","pvalue","padj")])

    
    ########################?????NEED to add limma ################################
    ########################?????NEED to add limma ################################
    #round(res,4)
    for(i in 1:4){
      res[,i] = data.frame(round(res[,i],4))
    }
    # for(i in 5:6){
    #   res[,i] = formatC(res[,i], format = "e", digits = 2)
    # }
  res
  #  data.frame(round(res,4))
  }
  # }
  #
  # , error=function(e) {
  #   print('Error: The selected dataset does not exist')
  # })
  
}
