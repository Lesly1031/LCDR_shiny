

deg_figUI<-function(headings,output_info,download_info,fill){
  column(6,
  panel(
    status="info",heading=paste0("Study: ",headings),
    dropdown(
      selectizeInput(fill, 'Optionally, input interested genes to highlight',
                     choices =NULL,multiple=T),
      downloadButton(download_info, "Download"),
      circle = TRUE, status = "danger", icon = icon("gear"), width = "300px",
      tooltip = tooltipOptions(title = "Click to see inputs !")
    ),
    conditionalPanel( 
      condition="input.update_MA==0",
      br(),
      p("Please click 'Run' button to generate results")
    ),
    
    conditionalPanel( 
      condition="input.update_MA>0",
      br(),
      withLoader(plotlyOutput(output_info,width = "100%",height = "300px"),type="html",loader="loader5")
      
      
    )
    

  ))
  
}

deg_fig_observe = function(session,input_info,deg_res){
  deg_res = data.frame(deg_res)
    updateSelectizeInput(session,input_info,choices  = unique(rownames(deg_res)),server=T)
  
}


############ deg visualize server ################
deg_figServer<-function(input,output,session,diff_result,folder_direct,fill_ma,title){
  
  
  
  if(input[["plot_type_MA"]] =="MA plot"){
    #   tryCatch({
    
    res=diff_result
    for(i in 1:6){
      res[,i] = as.numeric(as.character(res[,i]))
    }

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
      res$rowID = rownames(res)
      rownames(res) = NULL
      #res = data.frame(res[,c("rowID",names(res)[-ncol(res)]  )])
      if(grepl("metabolomics",folder_direct)|grepl("lipidomics",folder_direct)){ ############################################# change on server
        # names = sapply(strsplit(res$rowID,"_"),"[",1) 
        # IDs = unlist(lapply(1:length(res$rowID),function(x){
        #   paste(strsplit(res$rowID[x],"_")[[1]][2:3],collapse = "_")
        # }))
        split_sample_names <- function(names) {
          # The first part up to the polarity (pos or neg)
          part1 <- str_extract(names, ".*(?=_pos|_neg)")
          # The second part starting from the polarity (pos or neg)
          part2 <- str_extract(names, "(pos|neg)_.+")
          # Combine the parts into a data frame
          data.frame(part1, part2)
        }
        
        # Apply the function to split the names
        split_names_df <- split_sample_names(res$rowID)
        names = split_names_df$part1
        IDs = split_names_df$part2
        
        res$Names= ifelse(names=="X",res$rowID,names)
        res$IDs= ifelse(names=="X",res$rowID,IDs)
      }else{
        res$Names = sapply(strsplit(res$rowID,"_"),"[",1)
        res$IDs = sapply(strsplit(res$rowID,"_"),"[",-1)
      }
      
      if(input$which_p=="pval"){
        res$group=if_else(as.numeric(as.character(res$pval))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))
                          ,"Up",if_else(as.numeric(as.character(res$pval))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)<(-log2(input$fc_threshold_MA)),"Down","NS")) 
      }else{
        res$group=if_else(as.numeric(as.character(res$padj))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)>log2(as.numeric(input$fc_threshold_MA))
                          ,"Up",if_else(as.numeric(as.character(res$padj))<=as.numeric(input$p_threshold_MA)&as.numeric(res$log2FoldChange)<(-log2(input$fc_threshold_MA)),"Down","NS")) 
      }
      res$sel=ifelse(res$rowID%in%as.character(input[[fill_ma]]),"Selected","Not selected")
      res$group=factor(res$group,levels=c("Up","Down","NS"))
      
      res$log2_baseMean=log2(res$baseMean)
      res$gene_name <- res$Names
      res$delabel = NA
      res$delabel[res$sel=="Selected"] = res$Names[res$sel=="Selected"]
      length_res_delebel =sum(res$sel=="Selected")
      
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
    #   # load("F:/Projects/RNA_seq_shiny/RNAseq_shiny/rshiny/RNA_seq_shiny/temp_diff.RData")
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
      for(i in 1:6){
        res[,i] = as.numeric(as.character(res[,i]))
      }
      
      ###Figure begin###
      # load("F:/Projects/RNA_seq_shiny/RNAseq_shiny/rshiny/RNA_seq_shiny/temp_diff.RData")
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
        res$rowID = rownames(res)
        rownames(res) = NULL
        #res = data.frame(res[,c("rowID",names(res)[-ncol(res)]  )])
        if(grepl("metabolomics",folder_direct)|grepl("lipidomics",folder_direct)){
          # names = sapply(strsplit(res$rowID,"_"),"[",1) 
          # IDs = unlist(lapply(1:length(res$rowID),function(x){
          #   paste(strsplit(res$rowID[x],"_")[[1]][2:3],collapse = "_")
          # }))
          # res$Names= ifelse(names=="X",res$rowID,names)
          # res$IDs= ifelse(names=="X",res$rowID,IDs)
          
          split_sample_names <- function(names) {
            # The first part up to the polarity (pos or neg)
            part1 <- str_extract(names, ".*(?=_pos|_neg)")
            # The second part starting from the polarity (pos or neg)
            part2 <- str_extract(names, "(pos|neg)_.+")
            # Combine the parts into a data frame
            data.frame(part1, part2)
          }
          
          # Apply the function to split the names
          split_names_df <- split_sample_names(res$rowID)
          names = split_names_df$part1
          IDs = split_names_df$part2
          
          res$Names= ifelse(names=="X",res$rowID,names)
          res$IDs= ifelse(names=="X",res$rowID,IDs)
          
        }else{
          res$Names = sapply(strsplit(res$rowID,"_"),"[",1)
          res$IDs = sapply(strsplit(res$rowID,"_"),"[",-1)
          
        }
        
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
        res$gene_name <- as.factor(res$Names)
        res$sel=ifelse(res$rowID%in%c(input[[fill_ma]]),"Selected","Not selected")
        res$delabel = NA
        res$delabel[res$sel=="Selected"] = res$Names[res$sel=="Selected"]
        length_res_delebel =sum(res$sel=="Selected")
        #res$delabel[1:length_res_delebel] = sapply( strsplit( res$delabel[1:length_res_delebel], "_"), "[", 1) 

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
