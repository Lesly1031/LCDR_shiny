##############Heatmap##############
heatmapDegUI<-function(headings,output_info,color_annot){
  column(10,
         panel(
           status="info",heading=paste0("Study: ",headings),
           
           conditionalPanel( 
             condition="input.update_heat_deg==0",
             br(),
             p("Please click 'Run' button to generate figures, DEG tab has to be run prior to this tab")
           ),
           
           conditionalPanel( 
             condition="input.update_heat_deg>0",
             br(),
             textOutput(color_annot),
             (div(style='max-width: 100%; height: 500px; width: 90%;overflow-y: scroll;overflow-x: scroll;',
                  div(style="text-align: center;",plotOutput(output_info,width="100%",height="500px"))
                  
             ))              
           )
           
         )
  )
  
  
}



heatmapDeg_server = function(input,output,session,data,data_class,group_info,deg_res,group_info_deg){
  
    dt=data.frame(data,check.names = F)
    
  
  rownames(deg_res) = deg_res$gene_name
  #dt[is.na(dt)]<-0
  
  dt_class = data.frame(data_class)
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
  
  valid_columns <- dt_class$X %in% names(dt)
  selected_columns <- dt_class$X[valid_columns]
  
  dt <- dt %>% dplyr::select("ID",all_of(selected_columns))
  #dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names = F)
  
  
  dt_sub_genetype=dt
  sel_genes = rownames(deg_res)[1:input$top_n_heat_deg]
  dt_sub=subset(dt_sub_genetype,dt_sub_genetype$ID%in%sel_genes)
  
  dt_sub = data.frame(dt_sub,check.names = F)
  
  rownames(dt_sub)=dt_sub$ID
  dt_sub=data.frame(t(dt_sub[,-c(1)]),check.names = F)
    
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
    
  
  cal_z_score <- function(x){
    (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
  }
  dt_sub=apply(dt_sub, 2, cal_z_score)
  
  dt_class = data.frame(data_class)
  dt_class = dt_class[match(rownames(dt_sub),dt_class$X),]
  dt_class = dt_class[rowSums(is.na(dt_class)) != ncol(dt_class), ]
  
  group_info1 = str_split(input[[group_info]],":")[[1]][1]
  comp_group = str_split(input[[group_info]],":")[[1]][2]
  comp_group_sub = str_split(comp_group," vs ")[[1]]
  
  annotation = data.frame(dt_class[,group_info1])
  rownames(annotation) = dt_class$X
  names(annotation) = "group"
  annotation = subset(annotation,annotation$group%in%comp_group_sub)
  dt_sub =dt_sub[rownames(dt_sub)%in%rownames(annotation),]
  names_id = colnames(dt_sub)
  dt_sub=data.frame(t(dt_sub),check.names = F)
  dt_sub$status = status
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
  if(ncol(dt_sub)==1){
    
    if(input$show_row_names_deg){
      p= pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                  fontsize=font_size,show_rownames = T, show_colnames =T,na_col = "gray")
      
      cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
      #p$gtable$grobs[[5]]$gp=gpar(col=cols)
      p
    }else{
      pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               fontsize=font_size,show_rownames = input$show_row_names_deg, show_colnames = T,na_col = "gray")
      
      
      
    }
    
    # pheatmap(dt_sub,annotation = annotation,cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
    #          show_rownames =input$show_row_names_deg, show_colnames = input$show_column_names_deg,main="z-score heatmap")
  }else{
    if(input$show_row_names_deg){
      p= pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                  fontsize=font_size,show_rownames = T, show_colnames = T,na_col = "gray")
      
      cols=dt_sub[order(match(rownames(dt_sub), p$gtable$grobs[[5]]$label)), ]$colors
      #p$gtable$grobs[[5]]$gp=gpar(col=cols)
      p
    }else{
      pheatmap(dt_sub[,1:(ncol(dt_sub)-2)], annotation = annotation,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               fontsize=font_size,show_rownames = input$show_row_names_deg,show_colnames = T,
               main="z-score heatmap",na_col = "gray")
      
    }
    
    
    # pheatmap(dt_sub,annotation = annotation,cluster_rows=T,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
    #          show_rownames =input$show_row_names_deg, show_colnames = input$show_column_names_deg,main="z-score heatmap")
    
  }
  
}
