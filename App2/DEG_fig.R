############ deg visualize server ################
deg_figServer<-function(input,output,session,diff_result,folder_direct,fill_ma,title){
  
  
  
  if(input[["plot_type_MA"]] =="MA plot"){
    #   tryCatch({
    
    res=diff_result
    
    for(i in 1:ncol(res)){
      res[,i] = as.numeric(as.character(res[,i]))
    }
    names(res)[5]="pval"
    if(is.null(input[[fill_ma]])){
      if(input$which_p=="pval"){
        res$group=if_else(as.numeric(as.character(res$pval))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))
                          ,"Up",if_else(as.numeric(as.character(res$pval))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)<(-log2(input$fc_threshold_MA)),"Down","NS")) 
      }else{
        res$group=if_else(as.numeric(as.character(res$padj))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))
                          ,"Up",if_else(as.numeric(as.character(res$padj))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)<(-log2(input$fc_threshold_MA)),"Down","NS")) 
      }
      
      res$group=factor(res$group,levels=c("Up","Down","NS"))
      
      res$log2_baseMean=log2(res$baseMean)
      res$gene_name <- rownames(res)
      p=ggplot(res,aes_string(x="log2_baseMean",y="log2FoldChange",colour= "group",
                              text ="gene_name"))+
        geom_point()+
        labs(color = "")+
        scale_color_manual(values=c("coral2","#56B4E9","#999999"))+
        geom_hline(yintercept=0,linetype="dashed", color = "black")+
        theme_minimal()       +
        ggtitle(title)
    } else{
      if(input$which_p=="pval"){
        res$group=if_else(as.numeric(as.character(res$pval))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))
                          ,"Up",if_else(as.numeric(as.character(res$pval))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)<(-log2(input$fc_threshold_MA)),"Down","NS")) 
      }else{
        res$group=if_else(as.numeric(as.character(res$padj))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))
                          ,"Up",if_else(as.numeric(as.character(res$padj))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)<(-log2(input$fc_threshold_MA)),"Down","NS")) 
      }
      res$sel=ifelse(rownames(res)%in%as.character(make.names(input[[fill_ma]])),"Selected","Not selected")
      res$group=factor(res$group,levels=c("Up","Down","NS"))
      
      res$log2_baseMean=log2(res$baseMean)
      res$gene_name <- rownames(res)
      res$delabel = NA
      res$delabel[res$sel=="Selected"] = rownames(res)[res$sel=="Selected"]
      length_res_delebel =sum(res$sel=="Selected")
      # res$delabel = sapply(1:length(res$sel),function(x){
      #   strsplit( res$delabel[x], "_")[[1]][1]
      # })
      # sapply( strsplit( res$delabel[1:length_res_delebel], "_"), "[", 1)      
      p=ggplot(res,aes_string(x="log2_baseMean",y="log2FoldChange",colour="sel",text="gene_name"))+
        geom_point()+
        labs(color = "")+
        geom_hline(yintercept=0,linetype="dashed", color = "black")+
        theme_minimal()+
        geom_text(
          label = res$delabel
        )+
        scale_color_manual(values=c("#D3D3D3","#FF0000"))+
        ggtitle(title)
      p
    }
    
    # }
    # load("F:/Projects/RNA_seq_shiny/RNAseq_shiny/rshiny/temp_diff.RData")
    
    
    
    
    
    
    #   p
    
    #  }
    # ,error=function(e){
    #    p=ggplot() + theme_void()+ggtitle("'Differential Gene Expression Analysis' tab is required to process prior to this tab")
    #  #  p
    #  })
  } else if (input[["plot_type_MA"]]=="Volcano plot"){
    
    # if(grepl("metabolomics",folder_direct)){
    #   #  tryCatch({
    #   res=diff_result
    #   
    #   ###Figure begin###
    #   # load("/temp_diff.RData")
    #  # names(res)[which(names(res)=="fc")]="log2FoldChange"
    #   res$log2FoldChange = log2(res$fc)
    #   res$log_p.adj=-10*log10(res$padj)
    #   thres_p=-10*log10(as.numeric(input$p_threshold_MA))
    #   thres_logfc=as.numeric(input$fc_threshold_MA)
    #   res$group=if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>as.numeric(input$fc_threshold_MA)|res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)<(-as.numeric(input$fc_threshold_MA))
    #                     ,"p and log2 fold", if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>=(-as.numeric(input$fc_threshold_MA))&as.numeric(res$log2FoldChange)<=as.numeric(input$fc_threshold_MA),"p",
    #                                                 if_else( res$log_p.adj>thres_p&as.numeric(res$log2FoldChange)>as.numeric(input$fc_threshold_MA)|as.numeric(res$log2FoldChange)<(-as.numeric(input$fc_threshold_MA)),"log2FC","NC")))
    #   res$group<-factor( res$group,levels=c("p and log2 fold","p","log2FC","NC"))
    #   res$gene_name <- as.factor(rownames(res))
    #   
    #   
    #   p=ggplot(res,aes_string(x="log2FoldChange",y="log_p.adj",colour="group",text="gene_name"))+
    #     geom_point()+
    #     ylab("-10Log10(p.adj)")+
    #     xlab("FoldChange")+
    #     geom_hline(yintercept=thres_p,linetype="dotted")+
    #     geom_vline(xintercept=thres_logfc,linetype="dotted")+
    #     geom_vline(xintercept=-thres_logfc,linetype="dotted")+
    #     scale_color_manual(values=c("dark green","gray","coral2","#56B4E9"))+
    #     theme_minimal()+
    #     theme(legend.position = "none")
    #   
    #   
    #   #}
    #   # , error=function(e){
    #   #   p=ggplot() + theme_void()+ggtitle("'Differential Gene Expression Analysis' tab is required to process prior to this tab")
    #   #   p
    #   # }
    #   #)
    # }else{
    #  tryCatch({
    res=diff_result
    for(i in 1:ncol(res)){
      res[,i] = as.numeric(as.character(res[,i]))
    }
    names(res)[5]="pval"
    
    ###Figure begin###
    # load("/temp_diff.RData")
    if(is.null(input[[fill_ma]])){
      if(input$which_p=="pval"){
        # res$log_p.adj=-10*log10(res$pval)
        # thres_p=-10*log10(as.numeric(input$p_threshold_MA))
        res$log_p.adj=-1*log10(res$pval)
        thres_p=-1*log10(as.numeric(input$p_threshold_MA))
        thres_logfc=log2(as.numeric(input$fc_threshold_MA))
        res$group=if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA)))
                          ,"p and log2 fold", if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>=(-log2(as.numeric(input$fc_threshold_MA)))&as.numeric(res$log2FoldChange)<=log2(as.numeric(input$fc_threshold_MA)),"p",
                                                      if_else( res$log_p.adj>thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA))),"log2FC","NC")))
      }else{
        res$log_p.adj=-1*log10(res$padj)
        thres_p=-1*log10(as.numeric(input$p_threshold_MA))
        thres_logfc=log2(as.numeric(input$fc_threshold_MA))
        res$group=if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA)))
                          ,"p and log2 fold", if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>=(-log2(as.numeric(input$fc_threshold_MA)))&as.numeric(res$log2FoldChange)<=log2(as.numeric(input$fc_threshold_MA)),"p",
                                                      if_else( res$log_p.adj>thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA))),"log2FC","NC")))
      }
      
      res$group<-factor( res$group,levels=c("p and log2 fold","p","log2FC","NC"))
      res$gene_name <- as.factor(rownames(res))
      
      
      p=ggplot(res,aes_string(x="log2FoldChange",y="log_p.adj",colour="group",text="gene_name"))+
        geom_point()+
        ylab(paste0("-Log10_",input$which_p))+
        geom_hline(yintercept=thres_p,linetype="dotted")+
        geom_vline(xintercept=thres_logfc,linetype="dotted")+
        geom_vline(xintercept=-thres_logfc,linetype="dotted")+
        scale_color_manual(values=c("dark green","gray","coral2","#56B4E9"))+
        theme_minimal()+
        theme(legend.position = "none")+
        ggtitle(title)
      
      
      #}
      # , error=function(e){
      #   p=ggplot() + theme_void()+ggtitle("'Differential Gene Expression Analysis' tab is required to process prior to this tab")
      #   p
      # }
      #)
      # }
    }else{
      if(input$which_p=="pval"){
        res$log_p.adj=-10*log10(res$pval)
        thres_p=-10*log10(as.numeric(input$p_threshold_MA))
        thres_logfc=log2(as.numeric(input$fc_threshold_MA))
        res$group=if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA)))
                          ,"p and log2 fold", if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>=(-log2(as.numeric(input$fc_threshold_MA)))&as.numeric(res$log2FoldChange)<=log2(as.numeric(input$fc_threshold_MA)),"p",
                                                      if_else( res$log_p.adj>thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA))),"log2FC","NC")))
      }else{
        res$log_p.adj=-10*log10(res$padj)
        thres_p=-10*log10(as.numeric(input$p_threshold_MA))
        thres_logfc=log2(as.numeric(input$fc_threshold_MA))
        res$group=if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA)))
                          ,"p and log2 fold", if_else(res$log_p.adj<=thres_p&as.numeric(res$log2FoldChange)>=(-log2(as.numeric(input$fc_threshold_MA)))&as.numeric(res$log2FoldChange)<=log2(as.numeric(input$fc_threshold_MA)),"p",
                                                      if_else( res$log_p.adj>thres_p&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))|as.numeric(res$log2FoldChange)<(-log2(as.numeric(input$fc_threshold_MA))),"log2FC","NC")))
      }
      
      res$group<-factor( res$group,levels=c("p and log2 fold","p","log2FC","NC"))
      res$gene_name <- as.factor(rownames(res))
      res$sel=ifelse(rownames(res)%in%as.character(make.names(input[[fill_ma]])),"Selected","Not selected")
      res$delabel = NA
      res$delabel[res$sel=="Selected"] = rownames(res)[res$sel=="Selected"]
      length_res_delebel =sum(res$sel=="Selected")
      #res$delabel[1:length_res_delebel] = sapply( strsplit( res$delabel[1:length_res_delebel], "_"), "[", 1) 
      # res$delabel = sapply(1:length(res$sel),function(x){
      #   strsplit( res$delabel[x], "_")[[1]][1]
      # })
      p=ggplot(res,aes_string(x="log2FoldChange",y="log_p.adj",colour="sel",text="gene_name"))+
        geom_point()+
        ylab(paste0("-Log10_",input$which_p))+
        geom_hline(yintercept=thres_p,linetype="dotted")+
        geom_vline(xintercept=thres_logfc,linetype="dotted")+
        geom_vline(xintercept=-thres_logfc,linetype="dotted")+
        scale_color_manual(values=c("dark green","gray","coral2","#56B4E9"))+
        theme_minimal()+
        theme(legend.position = "none")+
        geom_text(
          label = res$delabel
        )+
        scale_color_manual(values=c("#D3D3D3","#FF0000"))+
        ggtitle(title)
    }
    
    
    
  }
  # ggplotly(p, source = "select", tooltip = c("key"))
  p
}
