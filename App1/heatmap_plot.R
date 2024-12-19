##############Heatmap##############
heat_choice = function(phenotype){
  if(is.data.frame(phenotype) && length(phenotype) > 0) {
    phenotype = data.frame(phenotype)
    which_conti = sapply(names(phenotype),function(x) !is.numeric(phenotype[[x]]) )
    phenotype = data.frame(phenotype[names(which(which_conti))])
    names(phenotype)[-1]
  } else {
    return(character(0))  # Returns an empty character vector if the input is not as expected
  }
  
}

heatmapUI<-function(headings,output_info,input_info,group_info,phenotype,color_annot){
  column(10,
         panel(
           status="info",heading=paste0("Study: ",headings),
           dropdown(
             # select2Input(input_info, label="Choose a gene",
             #             choices =geneID),
             # pickerInput(input_info, label="Choose a gene", 
             #                     choices =geneID,multiple=T),
             # fileInput(input_info_file, 
             #           span(
             #             list(HTML("<p><abbr title='The uploaded data has to be one column csv file with geneNames on each row'>Upload your genelist (csv file)...</abbr></p>"))
             #           ),
             #           accept = c(
             #             '.csv'
             #           ))
             selectizeInput(input_info, label="Optionally, choose genes from the dropdown list",
                         choices =NULL,multiple=T),
             # selectizeInput(input_info, label="Choose a gene", 
             #                choices =NULL,multiple=T),
             selectizeInput(group_info, label="Select group", 
                            choices =c("All",heat_choice(phenotype))),
             #  downloadButton(download_info, "Download"),
             circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
             tooltip = tooltipOptions(title = "Click to see inputs !")
           ),
           
           conditionalPanel( 
             condition="input.update_heat==0",
             br(),
             p("Please click 'Run' button to generate figures")
           ),
           
           conditionalPanel( 
             condition="input.update_heat>0",
             br(),
             textOutput(color_annot),
             (div(style='max-width: 100%; height: 500px; width: 90%;overflow-y: scroll;overflow-x: scroll;',
                  div(style="text-align: center;",plotOutput(output_info,width="100%",height="500px"))
                  
             ))              
           )
           
           
           
         )
  )
  
}

# exclude rows with variance = 0
heat_data_func = function(input_data){
  input_data = input_data[which(apply(input_data[,-1], 1, var,na.rm=T) != 0),]
  input_data$ID
}


heatmap_gene_sel = function(session,input_info,genotype,genotype_class,study_type){
  
  observe({
    if(is.null(genotype) || is.null(genotype_class)){
      return()  # Exit if data is not available
    }
    if(grepl("lipidomics",study_type)){
      updateSelectizeInput(session=session, input_info, choices = list(
        "Collapse by class" = heat_data_func(genotype_class), 
        "Individual lipidomics" = heat_data_func(genotype)
      ),server=T)
    }else{
      updateSelectizeInput(session=session, input_info, choices = list("select omics"= heat_data_func(genotype)),server=T)
    }
    
  })
}

heat_color_annot = function(input,output,session,study_type){
  if(grepl("metabolomics",study_type)&input$show_row_names){
    paste("Red character metabolites are those with positive modes, and blue are metabolites with negative modes")
  }else{
    paste("")
  }
}

check_pos_neg <- function(name) {
  if (grepl("_pos_", name)) {
    return("pos")
  } else if (grepl("_neg_", name)) {
    return("neg")
  } else {
    return(NA)
  }
}


####################### Use pheatmap ################################################
heatmapserver<-function(input,output,session,data,data_class,input_gene,group_info,deg_res,folder_direct,data_lip,study_type){
  
  
  tryCatch({
    
    if(grepl("lipidomics",study_type)){
      common_cols <- intersect(names(data), names(data_lip))
      dt <- data.frame(rbind(data[, common_cols], data_lip[, common_cols]),check.names=F)
      
    }else{
      dt=data.frame(data,check.names = F)
      
    }
    # if(input$heat_meta=="Yes"){
    #   dt=data.frame(data_identify)
    #   
    # }else{
    #   dt=data.frame(data)
    #   
    # }
    
    dt_class = data.frame(data_class)
    
    
    dt_class$X[dt_class$X==""]=NA
    dt_class = dt_class[complete.cases(dt_class$X),]
    dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names = F)
    
    dt_sub_genetype=dt      
    
    #dt_sub=subset(dt_sub_genetype,dt_sub_genetype$ID%in%c("ENSMUSG00000006390.15_Elovl1","ENSMUSG00000006392.16_Med8","ENSMUSG00000006418.17_Rnf114"))
    if(!is.null(input[[input_gene]])){
      dt_sub=subset(dt_sub_genetype,dt_sub_genetype$ID%in%c(input[[input_gene]]))
      
    }else{
      sel_genes = rownames(deg_res)[1:25]
      dt_sub=subset(dt_sub_genetype,dt_sub_genetype$ID%in%sel_genes)
      
    }
    dt_sub = data.frame(dt_sub,check.names=F)
    
    rownames(dt_sub)=dt_sub$ID
    dt_sub=data.frame(t(dt_sub[,-c(1)]),check.names=F)
    
    if(grepl("metabolomics",folder_direct)|grepl("lipidomics",folder_direct)){
      status <- sapply(names(dt_sub), check_pos_neg)
      
      extract_name <- function(name) {
        sub("(_pos_|_neg_).*", "", name)
      }
      
      # Apply the function to each sample name to get the list of names
      extracted_names <- sapply(names(dt_sub), extract_name)
      
      names(dt_sub) = sapply(1:ncol(dt_sub),function(x){
        res = extracted_names[x]
        if(res=="X"){
          names_dt_sub = names(dt_sub)[x]
        } else{
          names_dt_sub = res
        }
      })
      
    }else{
      status = sapply(1:ncol(dt_sub),function(x){
        strsplit(names(dt_sub)[x],"_")[[1]][2]
      })
      
      names(dt_sub) = sapply(1:ncol(dt_sub),function(x){
        res = strsplit(names(dt_sub)[x],"_")[[1]][1]
        if(res=="X"){
          names_dt_sub = names(dt_sub)[x]
        } else{
          names_dt_sub = res
        }
      })
      
    }
    
    check_pos_neg <- function(name) {
      if (grepl("_pos_", name)) {
        return("pos")
      } else if (grepl("_neg_", name)) {
        return("neg")
      } else {
        return(NA)
      }
    }
    
    
    cal_z_score <- function(x){
      (x - mean(x,na.rm=T)) / sd(x,na.rm = T)
    }
    dt_sub=apply(dt_sub, 2, cal_z_score)
    
    if(input[[group_info]]=="All"){
      names_id = colnames(dt_sub)
      dt_sub = data.frame(t(dt_sub),check.names = F)
      dt_sub$status = status
      dt_sub$colors=ifelse(dt_sub$status=="pos","red","blue")
      font_size=ifelse(nrow(dt_sub)<=20,12,ifelse(nrow(dt_sub)>20&nrow(dt_sub)<=50,10,
                                                  ifelse(nrow(dt_sub)>50&nrow(dt_sub)<=80,9,
                                                         ifelse(nrow(dt_sub)>80&nrow(dt_sub)<=100,8,
                                                                ifelse(nrow(dt_sub)>10&nrow(dt_sub)<=130,7,5)))))
      #rownames(dt_sub) = ifelse(duplicated(names_id),paste0(names_id,"_",status),names_id)
      rownames_dt_sub = ifelse(duplicated(names_id),paste0(names_id,"_",status),names_id)
      rownames(dt_sub) =  make.unique(rownames_dt_sub)
      
      # length_nchar =  sapply(1:nrow(dt_sub),function(x){
      #   nchar(rownames(dt_sub)[x])
      # })
      # 
      # length_nchar_l = which(length_nchar>20)
      # rownames(dt_sub)[length_nchar_l] = make.names(paste0(substr(rownames(dt_sub)[length_nchar_l],1,20),"..."),unique=TRUE)
      # 
      # length_nchar2 =  sapply(1:ncol(dt_sub),function(x){
      #   nchar(colnames(dt_sub)[x])
      # })
      # 
      # length_nchar_l2 = which(length_nchar2>20)
      # colnames(dt_sub)[length_nchar_l2] = make.names(paste0(substr(colnames(dt_sub)[length_nchar_l2],1,20),"..."),unique=TRUE)
      
      if(nrow(dt_sub)==1){
        
        pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                 show_rownames =input$show_row_names, show_colnames = T,main="z-score heatmap",na_col = "gray")
      }else{
        p= pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],cluster_rows=T,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                    show_rownames =input$show_row_names, show_colnames = T,main="z-score heatmap",na_col = "gray")
        
        cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
        p$gtable$grobs[[6]]$gp=gpar(col=cols)
        p
        
      }
      
    }else{
      dt_class = data.frame(data_class,check.names=F)
      dt_class = dt_class[match(rownames(dt_sub),dt_class$X),]
      dt_class = dt_class[rowSums(is.na(dt_class)) != ncol(dt_class), ]
      
      annotation = data.frame(dt_class[,input[[group_info]]])
      rownames(annotation) = dt_class$X
      names(annotation) = "group"
      
      dt_sub =dt_sub[rownames(dt_sub)%in%rownames(annotation),]
      names_id = colnames(dt_sub)
      dt_sub=data.frame(t(dt_sub),check.names = F)
      dt_sub$status = status
      dt_sub$colors=ifelse(dt_sub$status=="pos","red","blue")
      # rownames(dt_sub) = ifelse(duplicated(names_id),paste0(names_id,"_",status),names_id)
      rownames_dt_sub = ifelse(duplicated(names_id),paste0(names_id,"_",status),names_id)
      rownames(dt_sub) =  make.unique(rownames_dt_sub)
      font_size=ifelse(nrow(dt_sub)<=20,12,ifelse(nrow(dt_sub)>20&nrow(dt_sub)<=50,10,
                                                  ifelse(nrow(dt_sub)>50&nrow(dt_sub)<=80,9,
                                                         ifelse(nrow(dt_sub)>80&nrow(dt_sub)<=100,8,
                                                                ifelse(nrow(dt_sub)>10&nrow(dt_sub)<=130,7,5)))))
      # length_nchar =  sapply(1:nrow(dt_sub),function(x){
      #   nchar(rownames(dt_sub)[x])
      # })
      # 
      # length_nchar_l = which(length_nchar>20)
      # rownames(dt_sub)[length_nchar_l] = make.names(paste0(substr(rownames(dt_sub)[length_nchar_l],1,20),"..."),unique=TRUE)
      # 
      # 
      # length_nchar2 =  sapply(1:ncol(dt_sub),function(x){
      #   nchar(colnames(dt_sub)[x])
      # })
      # 
      # length_nchar_l2 = which(length_nchar2>20)
      # colnames(dt_sub)[length_nchar_l2] =make.names(paste0(substr(colnames(dt_sub)[length_nchar_l2],1,20),"..."),unique=TRUE)
      
      
      if(nrow(dt_sub)==1){
        
        pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],annotation = annotation,cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                 show_rownames =input$show_row_names, show_colnames = T,main="z-score heatmap",na_col = "gray")
      }else{
        if(input$show_row_names){
          p= pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                      fontsize=font_size,show_rownames = T, show_colnames = T,na_col = "gray")
          
          cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
          p$gtable$grobs[[5]]$gp=gpar(col=cols)
          p
        }else{
          pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   fontsize=font_size,show_rownames = F,show_colnames = T,
                   main="z-score heatmap",na_col = "gray")
          
        }
        
        # p =  pheatmap(dt_sub[,1:(ncol(dt_sub)-2)],annotation = annotation,cluster_rows=T,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
        #          show_rownames =input$show_row_names, show_colnames = T,main="z-score heatmap",na_col = "gray")
        # cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
        # p$gtable$grobs[[6]]$gp=gpar(col=cols)
        # p
      }
    }
    
  }
  # , error=function(e){
  #   plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  #   text(x = 0.5, y = 0.5, paste0("Please run DEG or select genes"),
  #        cex = 1.6, col = "black")
  # }
  )
}

######################## Use complex heatmap #########################################
# heatmapserver<-function(input,output,session,data,data_class,input_gene,group_info){
# 
#   
#   tryCatch({
# 
#     dt=data.frame(data)
#     # for(i in 2:ncol(dt)){
#     #   dt[,i] = log2(dt[,i]+1)
#     # }
# 
#     dt_class=data.frame(data_class)
# 
# 
#     dt[is.na(dt)]<-0
# 
#     rownames_dt_class=as.character(dt_class[,1])
#     colnames_dt_class=colnames(dt_class)[-1]
#     dt_class=data.frame(dt_class[,-1])
#     rownames(dt_class)=rownames_dt_class
#     names(dt_class)=colnames_dt_class
# 
#     dt_sub_genetype=dt
# 
#     #dt_sub=subset(dt_sub_genetype,dt_sub_genetype$ID%in%c("ENSMUSG00000006390.15_Elovl1","ENSMUSG00000006392.16_Med8","ENSMUSG00000006418.17_Rnf114"))
#     
#     dt_sub=subset(dt_sub_genetype,dt_sub_genetype$ID%in%as.character(input[[input_gene]]))
#     dt_sub = data.frame(dt_sub)
#     
#     rownames(dt_sub)=dt_sub$ID
#     dt_sub=data.frame(t(dt_sub[,-c(1)]))
# 
# 
# 
#     # dt_sub=subset(dt_sub_genetype,dt_sub_genetype$ID%in%c("TSPAN6","C1orf112","FUCA2","GCLC","DPM1"))
# 
#       dt_class=data.frame(dt_class)
#       rownames(dt_class)=make.names(rownames(dt_class))
#       dt_sub=merge(dt_class,dt_sub,by="row.names")
#       rownames(dt_sub)=dt_sub$Row.names
#       dt_sub1=data.frame(dt_sub[,-1])
# 
#       if(input$z_score==TRUE){
#         cal_z_score <- function(x){
#           (x - mean(x)) / sd(x)
#         }
#         dt_sub=apply(dt_sub1[,-c(1:ncol(dt_class))], 2, cal_z_score)
#         dt_sub=merge(dt_class,dt_sub,by="row.names")
#         rownames(dt_sub)=dt_sub$Row.names
#         dt_sub=data.frame(dt_sub[,-1])
#       }else{
#         dt_sub=dt_sub1
#       }
# 
# 
# 
# 
#     dt_sub=data.frame(t(dt_sub))
#     dt_sub=data.frame(na.omit(dt_sub))
#     height_map=ifelse(nrow(dt_sub)<=20,8,ifelse(nrow(dt_sub)>20&nrow(dt_sub)<=50,10,
#                                                 ifelse(nrow(dt_sub)>50&nrow(dt_sub)<=80,15,
#                                                        ifelse(nrow(dt_sub)>80&nrow(dt_sub)<=100,18,
#                                                               ifelse(nrow(dt_sub)>10&nrow(dt_sub)<=130,22,30)))))
#     font_size=ifelse(nrow(dt_sub)<=20,12,ifelse(nrow(dt_sub)>20&nrow(dt_sub)<=50,10,
#                                                 ifelse(nrow(dt_sub)>50&nrow(dt_sub)<=80,9,
#                                                        ifelse(nrow(dt_sub)>80&nrow(dt_sub)<=100,8,
#                                                               ifelse(nrow(dt_sub)>10&nrow(dt_sub)<=130,7,5)))))
#     #names(dt_sub)=colnames(dt)[-1]
# 
# 
#       dt_sub=dt_sub[-c(1:ncol(dt_class)),]
# 
# 
#     for(i in 1: ncol(dt_sub)){
#       dt_sub[,i]=as.numeric(as.character(dt_sub[,i]))
# 
#     }
# 
#     if(input[[group_info]]=="All"){
#       Heatmap(dt_sub,show_row_names=input$show_row_names,show_column_names=input$show_column_names,
#               cluster_rows=input$cluster_rows,cluster_columns=input$cluster_columns,heatmap_legend_param = list(title = "Value"),height = unit(height_map, "cm")
#               ,row_names_gp=gpar(fontsize = font_size))
#     }else{
#       #names(dt_sub)=sample_name
#       #rownames(dt_class)=sample_name
#       dt_class = dt_class[match(names(dt_sub),rownames(dt_class)),]
#       # sample_anno=HeatmapAnnotation(group=as.character(dt_class[,which(names(dt_class)==input[[group_info]])]),
#       #                               col=degColors(as.character(dt_class[,which(names(dt_class)==input[[group_info]])]),FALSE,palette = "Dark2"))
#       sample_anno=HeatmapAnnotation(group=as.character(dt_class[,which(names(dt_class)==input[[group_info]])]))
#       Heatmap(dt_sub,show_row_names=input$show_row_names,show_column_names=input$show_column_names,
#               cluster_rows=input$cluster_rows,cluster_columns=input$cluster_columns,top_annotation = sample_anno,
#               heatmap_legend_param = list(title = "Value"), height = unit(height_map, "cm"),row_names_gp=gpar(fontsize = font_size))
#     }
#   }
#   , error=function(e){
#     plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#     text(x = 0.5, y = 0.5, paste0("Please select genes"),
#          cex = 1.6, col = "black")
#   }
#   )
# }
