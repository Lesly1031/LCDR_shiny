kegg_heatUI<-function(headings,output_info,output_info_tb,kegg_path,group_info,enrich_data,phenotype,download_info,gene_pathway,color_annot){
    column(10,
           panel(
             status="info",heading=paste0("Study: ",headings),
             dropdown(
               selectInput(kegg_path, h4("Select KEGG pathaway"), 
                           choices = NULL),
               selectInput(group_info, h4("Sample annotation"), 
                           choices = c(colnames(phenotype)[-1],"None"),selected = "None"),
               downloadButton(download_info, "Download"),
               circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
               tooltip = tooltipOptions(title = "Click to see inputs !")
             ),
             
             tabsetPanel(

               tabPanel("KEGG heatmap",
                        conditionalPanel( 
                          condition="input.update_heat_kegg>0",
                          br(),
                          textOutput(color_annot),
                          (div(style='max-width: 100%; height: 500px; width: 90%;overflow-y: scroll;overflow-x: scroll;',
                               div(style="text-align: center;",plotOutput(output_info,width="100%",height="500px"))
                               
                          ))              
                        )
                        ),
               tabPanel("KEGG table",
                        (div(style='max-width: 100%; height: 500px; width: 90%;overflow-y: scroll;overflow-x: scroll;',
                             div(style="text-align: center;",DTOutput(output_info_tb))
                             
                        ))   
               )

             ),

      
             conditionalPanel(
               condition="input.update_heat_kegg==0",
               br(),
               p("Select genes by clicking the gear button, then click 'Update' button to get the figure")
             )

             
           )) 
 
}

enrich_res_tb = function(input,output,session,deg_data,folder_direct){
 tryCatch({
   if(grepl("metabolomics",folder_direct)){
     res = deg_data
     res$gene_name <- rownames(res)
     for(i in 1:6){
       res[,i] = as.numeric(as.character(res[,i]))
     }
     res = data.frame(res[which(abs(res$pval)<0.05),])
     
     ID = sapply(1:length(res$gene_name ),function(x){
       str_split(res$gene_name[x], "_")[[1]][1]
     })
     ID = unique(ID)
    
     
     enriched = fread("/kegg_metabolites_processed.csv")
     enrich = data.frame(enriched)
     enrich = subset(enrich,enrich$Term!="")
     filter_genes <- function(gene_string, genes_to_keep) {
       split_genes <- str_split(gene_string, ";")[[1]]
       # Keep only the genes of interest
       filtered_genes <- split_genes[split_genes %in% genes_to_keep]
       # Combine them back into a single string
       paste(filtered_genes, collapse = ";")
     }
     
     # enrich$Genes = sapply(enrich$Genes, filter_genes, genes_to_keep = ID)
     # enrich = data.frame(enrich[order(enrich$length_genes,decreasing = T),])
     
   }else{
     res = deg_data
     res$gene_name <- rownames(res)
     for(i in 1:6){
       res[,i] = as.numeric(as.character(res[,i]))
     }
     res = data.frame(res[which(abs(res$pval)<0.05),])
     
     # if(input$which_p_kegg=="pval"){
     #   res = data.frame(res[which(abs(res$pval)<input$p_threshold_kegg&abs(res$log2FoldChange)>log2(input$fc_threshold_kegg)),])
     #   
     # }else if (input$which_p_kegg=="adjust.p"){
     #   res = data.frame(res[which(abs(res$padj)<input$p_threshold_kegg&abs(res$log2FoldChange)>log2(input$fc_threshold_kegg)),])
     #   
     # }
     
     # ID =  sapply(str_split(rownames_res, "_"), "[", 1)
     ID = sapply(1:length(res$gene_name ),function(x){
       str_split(res$gene_name[x], "_")[[1]][1]
     })
     ID = unique(ID)
     if(str_detect(ID[1], "^[:upper:]+$")){ # if all upper case - human data, else mice data
       dbs <- c("KEGG_2019_Human")
       enriched1 <- enrichr(ID, dbs)
       enrich_score = enriched1[["KEGG_2019_Human"]]
     }else{
       dbs <- c("KEGG_2019_Mouse")
       enriched1 <- enrichr(ID, dbs)
       enrich_score = enriched1[["KEGG_2019_Mouse"]]
     }
     enrich = data.frame(enrich_score[,c("Term",'Odds.Ratio',"Combined.Score","P.value","Adjusted.P.value","Genes")])
     # enrich = enrich_score
     enrich = data.frame(enrich[order(enrich$Adjusted.P.value,decreasing = F),])
   }
   
   enrich_genes =  strsplit(enrich$Genes, ";")
   enrich$length_genes = as.numeric(lengths(enrich_genes))
   enrich =subset(enrich,enrich$length_genes!=0)
   enrich <- enrich[ , !names(enrich) %in% "length_genes"]
   rownames(enrich) = NULL
   enrich
 },error = function(e){
   "No genes can be selected, please change the threshold"
   
 })


}

kegg_color_annot = function(input,output,session,folder_direct){
  if(grepl("metabolomics",folder_direct)&input$show_row_names_kegg){
    paste("Red character metabolites are those with positive modes, and blue are metabolites with negative modes")
  }else{
    paste("")
  }
}
kegg_sel_observe = function(session,input_info,enrich_res){
  tryCatch({

        updateSelectizeInput(session,input_info,choices = enrich_res$Term,selected = enrich_res$Term[5] ,server = T)
    
  },error = function(e){
    NULL
  })

}

#################################### pheatmap ##############################################################
kegg_heatServer<-function(input,output,session,enrich_data,omics_data,kegg_path,group_info,phenotype_data,folder_direct,group_info_deg){
  
 # tryCatch({
    enrich = data.frame(enrich_data)
    

    enrich_genes=enrich$Genes[enrich$Term==input[[kegg_path]]]
    enrich_genes <-  strsplit(enrich_genes, ";")[[1]]
    
    omics = data.frame(omics_data,check.names = F)
    dt_class = data.frame(phenotype_data)
    dt_class$X[dt_class$X==""]=NA
    dt_class = dt_class[complete.cases(dt_class$X),]
    
    
    # subset dt_class with only showing those that been searched in DEG
    input_string <- input[[group_info_deg]]
    split_parts <- str_split(input_string, ":| vs ")
    group_col <- split_parts[[1]][1]
    comparison_1 <- split_parts[[1]][2]
    comparison_2 <- split_parts[[1]][3]
    
    dt_class <- dt_class %>% 
      filter(!!sym(group_col) %in% c(comparison_1, comparison_2))
    
    valid_columns <- dt_class$X %in% names(omics)
    selected_columns <- dt_class$X[valid_columns]
    
    omics <- omics %>% dplyr::select("ID",all_of(selected_columns))
    
    omics = data.frame(omics,check.names=F)
    omics$geneName =  sapply( strsplit(omics$ID, "_"), "[", 1)
    
    if(grepl("metabolomics",folder_direct)){
      omics_sub = data.frame(omics[omics$geneName%in%enrich_genes,],check.names = F)
      
    }else{
      omics_sub = data.frame(omics[toupper(omics$geneName)%in%enrich_genes,],check.names = F)
      
    }
    omics_sub = data.frame(omics_sub[,-ncol(omics_sub)],check.names = F)
    
    rownames(omics_sub) = omics_sub$ID
    omics_sub = data.frame(omics_sub[,-1],check.names = F)
    omics_sub= data.frame(t(omics_sub),check.names=F)
    

    for(i in 1:ncol(omics_sub)){
      omics_sub[,i] = ( omics_sub[,i] - mean( omics_sub[,i])) / sd( omics_sub[,i])
    }
    
    if(ncol(omics_sub)==1){
      dt_sub = omics_sub
    }else{
      dt_sub = omics_sub[ , which(apply(omics_sub, 2, var) != 0)]
      
    }
    
    height_map=ifelse(ncol(dt_sub)<=20,8,ifelse(ncol(dt_sub)>20&ncol(dt_sub)<=50,10,
                                                ifelse(ncol(dt_sub)>50&ncol(dt_sub)<=80,15,
                                                       ifelse(ncol(dt_sub)>80&ncol(dt_sub)<=100,18,
                                                              ifelse(ncol(dt_sub)>10&ncol(dt_sub)<=130,22,30)))))
    font_size=ifelse(ncol(dt_sub)<=20,12,ifelse(ncol(dt_sub)>20&ncol(dt_sub)<=50,10,
                                                ifelse(ncol(dt_sub)>50&ncol(dt_sub)<=80,9,
                                                       ifelse(ncol(dt_sub)>80&ncol(dt_sub)<=100,8,
                                                              ifelse(ncol(dt_sub)>10&ncol(dt_sub)<=130,7,5)))))
    
    status = sapply(1:ncol(dt_sub),function(x){
      strsplit(names(dt_sub)[x],"_")[[1]][2]
    })
     names(dt_sub) = sapply(1:ncol(dt_sub),function(x){
      strsplit(names(dt_sub)[x],"_")[[1]][1]
    })
    
    
    
    if(input[[group_info]]=="None"){
      if(ncol(dt_sub)==1){
        original_rownames = names(dt_sub) ##############
        original_rownames = ifelse(duplicated(original_rownames),paste0(original_rownames,"_",status),original_rownames) #############
        dt_sub = data.frame(t(dt_sub),check.names=F)
        dt_sub$status = status
        dt_sub$colors=ifelse(dt_sub$status=="pos","red","blue")
        length_nchar =  sapply(1:nrow(dt_sub),function(x){
          nchar(rownames(dt_sub)[x])
        })
        rownames(dt_sub) = original_rownames
        
        # length_nchar_l = which(length_nchar>10)
        # rownames(dt_sub)[length_nchar_l] = paste0(substr(rownames(dt_sub)[length_nchar_l],1,10),"...")
        # 
        # length_nchar2 =  sapply(1:ncol(dt_sub),function(x){
        #   nchar(colnames(dt_sub)[x])
        # })
        # 
        # length_nchar_l2 = which(length_nchar2>15)
        # colnames(dt_sub)[length_nchar_l2] = paste0(substr(colnames(dt_sub)[length_nchar_l2],1,15),"...")
        if(input$show_row_names_kegg){
          # p=pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
          #            show_rownames = T, show_colnames = T)
          # cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
          # p$gtable$grobs[[5]]$gp=gpar(col=cols)
          # p
          pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                   show_rownames = input$show_row_names_kegg, show_colnames = T)
        }else{
          pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                   show_rownames = input$show_row_names_kegg, show_colnames = T)
        }

      }else{
        original_rownames = names(dt_sub) ##############
        original_rownames = ifelse(duplicated(original_rownames),paste0(original_rownames,"_",status),original_rownames) #############
        
        dt_sub = data.frame(t(dt_sub),check.names=F)
        rownames(dt_sub) = original_rownames
        dt_sub$status = status
        dt_sub$colors=ifelse(dt_sub$status=="pos","red","blue")
        if(input$show_row_names_kegg){
          p=pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                     show_rownames = T, show_colnames = T)
          cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
          p$gtable$grobs[[5]]$gp=gpar(col=cols)
          p
        }else{
          pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                   show_rownames = input$show_row_names_kegg, show_colnames = T)
        }


      }

    }else{
      
      #########################
      dt_class = data.frame(phenotype_data,check.names=F)
      dt_class = dt_class[match(rownames(dt_sub),dt_class$X),]
      dt_class = dt_class[rowSums(is.na(dt_class)) != ncol(dt_class), ]
      #################################
      

      annotation = data.frame(dt_class[,input[[group_info]]],check.names=F)
      rownames(annotation) = dt_class$X
      names(annotation) = "group"
      
      original_rownames = names(dt_sub) ##############
      original_rownames = ifelse(duplicated(original_rownames),paste0(original_rownames,"_",status),original_rownames) #############
      
      dt_sub = data.frame(t(dt_sub),check.names=F)
      rownames(dt_sub) = original_rownames
      
      dt_sub$status = status
      dt_sub$colors=ifelse(dt_sub$status=="pos","red","blue")
      length_nchar =  sapply(1:nrow(dt_sub),function(x){
        nchar(rownames(dt_sub)[x])
      })

      # length_nchar_l = which(length_nchar>10)
      # rownames(dt_sub)[length_nchar_l] = paste0(substr(rownames(dt_sub)[length_nchar_l],1,10),"...")
      # 
      # length_nchar2 =  sapply(1:ncol(dt_sub),function(x){
      #   nchar(colnames(dt_sub)[x])
      # })
      # 
      # length_nchar_l2 = which(length_nchar2>15)
      # colnames(dt_sub)[length_nchar_l2] = paste0(substr(colnames(dt_sub)[length_nchar_l2],1,15),"...")

     if(ncol(dt_sub)==1){
       if(input$show_row_names_kegg){
         p= pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             fontsize=font_size,show_rownames = T, show_colnames =T,cluster_rows = FALSE)
         
         cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
         p$gtable$grobs[[5]]$gp=gpar(col=cols)
         p
       }else{
         pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                  fontsize=font_size,show_rownames = input$show_row_names_kegg, show_colnames = T,cluster_rows = FALSE)
         

         
       }

     } else{
       if(input$show_row_names_kegg){
         p= pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                            fontsize=font_size,show_rownames = T, show_colnames = T,cluster_rows=F)
         
         #cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
         #p$gtable$grobs[[5]]$gp=gpar(col=cols)
         p
       }else{
         pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                  fontsize=font_size,show_rownames = input$show_row_names_kegg, show_colnames = T,cluster_rows=F)
         
       }

     }

    }
  # },error = function(e){
  #   plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  #   text(x = 0.5, y = 0.5, paste0("No genes in the selected pathway"),
  #        cex = 1.6, col = "black")
  # })
  
  
}
################################ complexheatmap ##############################################################
# kegg_heatServer<-function(input,output,session,enrich_data,omics_data,kegg_path,group_info,phenotype_data){
# 
#   tryCatch({
#     enrich = data.frame(enrich_data)
#     enrich_genes=enrich$Genes[enrich$Term==input[[kegg_path]]]
#     enrich_genes <-  strsplit(enrich_genes, ";")[[1]]
#     
#     omics = data.frame(omics_data)
#     omics$geneName =  sapply( strsplit(omics$ID, "_"), "[", 1)
#     omics_sub = data.frame(omics[toupper(omics$geneName)%in%enrich_genes,])
#     omics_sub = data.frame(omics_sub[,-ncol(omics_sub)])
#     
#     rownames(omics_sub) = omics_sub$ID
#     omics_sub = data.frame(omics_sub[,-1])
#     omics_sub= data.frame(t(omics_sub))
#     
#     for(i in 1:ncol(omics_sub)){
#       omics_sub[,i] = ( omics_sub[,i] - mean( omics_sub[,i])) / sd( omics_sub[,i])
#     }
#     
#     dt_sub = omics_sub[ , which(apply(omics_sub, 2, var) != 0)]
#     height_map=ifelse(ncol(dt_sub)<=20,8,ifelse(ncol(dt_sub)>20&ncol(dt_sub)<=50,10,
#                                                 ifelse(ncol(dt_sub)>50&ncol(dt_sub)<=80,15,
#                                                        ifelse(ncol(dt_sub)>80&ncol(dt_sub)<=100,18,
#                                                               ifelse(ncol(dt_sub)>10&ncol(dt_sub)<=130,22,30)))))
#     font_size=ifelse(ncol(dt_sub)<=20,12,ifelse(ncol(dt_sub)>20&ncol(dt_sub)<=50,10,
#                                                 ifelse(ncol(dt_sub)>50&ncol(dt_sub)<=80,9,
#                                                        ifelse(ncol(dt_sub)>80&ncol(dt_sub)<=100,8,
#                                                               ifelse(ncol(dt_sub)>10&ncol(dt_sub)<=130,7,5)))))
#     
#     
#     if(input[[group_info]]=="None"){
#       Heatmap(as.matrix(t(dt_sub)),show_row_names=input$show_row_names_kegg,show_column_names=input$show_column_names_kegg,
#               cluster_rows=T,cluster_columns=T,heatmap_legend_param = list(title = "Value"),height = unit(8, "cm"),row_names_gp=gpar(fontsize = font_size))
#       
#     }else{
#       dt_class = data.frame(phenotype_data)
#       dt_class$X = make.names(dt_class$X)
#       dt_class = dt_class[match(rownames(dt_sub),dt_class$X),]
#       sample_anno=HeatmapAnnotation(group=as.character(dt_class[,which(names(dt_class)==input[[group_info]])]))
#       Heatmap(as.matrix(t(dt_sub)),show_row_names=input$show_row_names_kegg,show_column_names=input$show_column_names_kegg,
#               cluster_rows=T,cluster_columns=T,top_annotation = sample_anno,
#               heatmap_legend_param = list(title = "Value"), height = unit(8, "cm"),row_names_gp=gpar(fontsize = font_size))
#     }
#   },error = function(e){
#     plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#     text(x = 0.5, y = 0.5, paste0("No genes in the selected pathway"),
#          cex = 1.6, col = "black")
#   })
# 
# 
# }