pca_choice_shape = function(phenotype){
  # phenotype = data.frame(phenotype)
  # which_conti = sapply(names(phenotype),function(x) !is.numeric(phenotype[[x]]) )
  # phenotype = data.frame(phenotype[names(which(which_conti))])
  # names(phenotype)[-1]
  phenotype = data.frame(phenotype)
  
  # Identify variables that are not numeric and have 10 or fewer unique values
  which_valid = sapply(names(phenotype), function(x) {
    if (!is.numeric(phenotype[[x]])) {
      return(length(unique(phenotype[[x]])) <= 10)  # Check if unique values <= 10
    }
    return(FALSE)  # Return false for numeric (continuous) variables
  })
  
  # Select only those variables that meet the criteria (non-numeric and <= 10 unique values)
  phenotype = data.frame(phenotype[names(which(which_valid))])
  
  # Return the names of the selected variables
  names(phenotype)
  
}


pca_UI<-function(headings,output_info,group_info1,group_info2,download_info,phenotype,study_type){

    column(10,
           panel(
             status="info",heading=paste0("Study: ",headings),
             dropdownButton(
               selectizeInput(group_info1, label="Select color annotation", 
                              choices=c(names(phenotype)[-1])),   
               selectizeInput(group_info2, label="Select shape annotation", 
                              choices=c("None",pca_choice_shape(phenotype)[-1]),selected = "None"),    
               downloadButton(download_info, "Download"),
               circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
               tooltip = tooltipOptions(title = "Click to see inputs !")
             ),
             
             conditionalPanel( 
               condition="input.update_pca==0",
               br(),
               p("Please click 'Run' button to generate figures")
             ),
             
             conditionalPanel( 
               condition="input.update_pca>0",
               br(),
               withLoader(plotlyOutput(output_info),type="html",loader="loader5")
               #DTOutput(output_info)
             )
             
             
           )
           
    )
}

#####################################PCA functions###############################################################
##############################PCA functions###################################################
# data_class = read_fst("F:/Projects/Brooke/lung_repo/Shiny/shiny_explore030321/studies/5/metabolomics/phenotype_data.fst")
# dt = read_fst("F:/Projects/Brooke/lung_repo/Shiny/shiny_explore030321/studies/5/metabolomics/genotype_data.fst")

pca_server<-function(input,output,session,data,data_class,group_info1,group_info2,study_type){
  
  dt=data.frame(data,check.names = F)
  
  # for(i in 2:ncol(dt)){
  #   dt[,i] = log2(dt[,i]+1)
  # }
  
  dt_class=data.frame(data_class,check.names = F)
  
  dt_class$X[dt_class$X==""]=NA
  dt_class = dt_class[complete.cases(dt_class$X),]
  dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names = F)
  
  
  dt_sub=dt
  
  sample_name=colnames(dt_sub)[-1]
  dt_sub1=data.frame(t(dt_sub[,-c(1)]),check.names = F)
  rownames(dt_sub1)=sample_name
  colnames(dt_sub1)=dt_sub$ID
  dt_sub=dt_sub1
  #
  rownames_dt_class = as.character(na.omit(dt_class$X))
  rownames_dt_class[rownames_dt_class==""]=NA
  rownames_dt_class=na.omit(rownames_dt_class)
  colnames_dt_class=colnames(dt_class)[-1]
  dt_class=data.frame(dt_class[,-1],check.names = F)
  rownames(dt_class)=rownames_dt_class
  names(dt_class)=colnames_dt_class
  rownames(dt_class) = rownames(dt_class)
  
  #---------------- PCA ----------------#
  tryCatch({
    
    
    # if(input[[group_info1]]=="All"&input[[group_info2]]=="All"){
    #   
    #   
    #   
    #   dt_sub = dt_sub[ , which(apply(dt_sub, 2, var,na.rm=T) != 0)]
    #   
    #   if(grepl("metabolomics",folder_direct)|grepl("proteomics",folder_direct)|grepl("lipidomics",folder_direct)){
    #     nipals_result <- nipals(dt_sub, scale= TRUE, center=TRUE,ncomp = 2)
    #     scores = data.frame(nipals_result$scores,check.names = F)
    #     scores$X = rownames(scores)
    #     
    #     
    #     p= ggplot(scores,aes(x=PC1,y=PC2,text=X
    #     ))+
    #       geom_point()+
    #       theme_light()+
    #       xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
    #       ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))
    #     
    #   } else{
    #     
    #     pca_res <- prcomp(na.omit(dt_sub),center=T,scale.=T)
    #     
    #     t=data.frame(pca_res$x,check.names = F)
    #     scores = data.frame(X=rownames(t),PC1 = t$PC1, PC2 = t$PC2,check.names = F)
    #     
    #     p= ggplot(scores,aes(x=PC1,y=PC2,text=X
    #     ))+
    #       geom_point()+
    #       theme_light()+
    #       xlab(paste0("PC1"," (VarExp:",round(summary(pca_res)$importance[2,1],4)*100,"%)"))+
    #       ylab(paste0("PC2"," (VarExp:",round(summary(pca_res)$importance[2,2],4)*100,"%)"))
    #     
    #   }
    #   p
    #   #ggplotly(p, tooltip="text")
    # }
    # else if(input[[group_info]]=="eligible samples (For metabolomics)"){
    #         names_check <- sapply( strsplit( colnames(dt_sub), "_"), "[", 1)
    #         keep_index = which(names_check!="X")
    #         dt_sub = dt_sub[,keep_index]
    #
    #         nipals_result <- nipals(dt_sub, scale= TRUE, center=TRUE,ncomp = 2)
    #         scores = data.frame(nipals_result$scores)
    #         scores$X = rownames(scores)
    #
    #
    #         p= ggplot(scores,aes(x=PC1,y=PC2,text=X
    #         ))+
    #           geom_point()+
    #           theme_light()+
    #           xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
    #           ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))
    #         p
    #
    #     }
    #else{
      # # rownames_dt_class=as.character(dt_class[,1])
      # # colnames_dt_class=colnames(dt_class)[-1]
      # # dt_class=data.frame(dt_class[,-1])
      # # rownames(dt_class)=rownames_dt_class
      # # names(dt_class)=colnames_dt_class
      # dt_sub[is.na(dt_sub)]<-0
      #
      # dt_sub=merge(dt_class,dt_sub,by="row.names")
      # rownames(dt_sub)=dt_sub$Row.names
      # dt_sub=data.frame(dt_sub[,-1])
      # rownames(dt_class)=make.names(rownames(dt_class))
      # dt_sub[is.na(dt_sub)]<-0
      # pca_res <- prcomp(na.omit(dt_sub[,-c(1:ncol(dt_class))]))
      # t=data.frame(pca_res$x)
      # scores = data.frame(X=rownames(t),PC1 = t$PC1, PC2 = t$PC2,dt_sub[,1:ncol(dt_class)])
      #
      # colnames(scores)[4:ncol(scores)]=names(dt_class)
      # #  colnames(scores)[3:ncol(scores)]=names(dt_class)
      # for(i in 2:3){
      #   scores[,i]=as.numeric(as.character(scores[,i]))
      # }
      #
      # for(i in 4:ncol(scores)){
      #   scores[,i]=as.character(scores[,i])
      # }

      if(grepl("RNAseq",study_type)|grepl("Metabolomics",study_type)|grepl("Proteomics",study_type)|grepl("Lipidomics",study_type)){
        dt_sub = dt_sub[ , which(apply(dt_sub, 2, var,na.rm=T) != 0)]
        
        nipals_result <- nipals(dt_sub, scale= TRUE, center=TRUE,ncomp = 2)
        scores = data.frame(nipals_result$scores,check.names = F)
        scores$X = rownames(scores)
        scores = merge(scores,dt_class,by.x="X",by.y="row.names")
        for(i in 2:3){
          scores[,i]=as.numeric(as.character(scores[,i]))
        }
        
       #  if(input[[group_info1]]=="None"&input[[group_info2]]=="None"){
       #    p= ggplot(scores,aes(x=PC1,y=PC2,text = X))+
       #      geom_point(size = 2)+
       #      theme_light()+
       #      xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
       #      ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))
       #  }else if(input[[group_info1]]!="None"&input[[group_info2]]=="None"){ # if annotate by color not by shape
       #   if(is.numeric(scores[[input[[group_info1]]]])){ # if the color annotation is a continuous variable
       #     p= ggplot(scores,aes(x=PC1,y=PC2,text = X
       #     ))+
       #       geom_point(aes_string(colour = input[[group_info1]]),size = 2)+
       #       theme_light()+
       #       scale_color_gradient(low = "blue", high = "red") +
       #       xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
       #       ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))
       #   }else{ # if the color annotation is not a continuous variable
       #     p= ggplot(scores,aes(x=PC1,y=PC2,text = X
       #     ))+
       #       geom_point(aes_string(colour = input[[group_info1]]),size = 2)+
       #       theme_light()+
       #       xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
       #       ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))+
       #       labs(colour="")
       # 
       #   }
       #   
       # 
       # }else if(input[[group_info1]]=="None"&input[[group_info2]]!="None"){ # if annotate by shape not by color
       #   p= ggplot(scores,aes(x=PC1,y=PC2,text = X
       #   ))+
       #     geom_point(aes_string(shape = input[[group_info2]]),size = 2)+
       #     theme_light()+
       #     xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
       #     ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))+
       #     scale_shape_manual(values=1:200)+
       #     labs(shape="")
       #   
       # }
        if(input[[group_info2]]=="None"){ # if annotate by color not by shape
          if(is.numeric(scores[[input[[group_info1]]]])){ # if the color annotation is a continuous variable
            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(colour = input[[group_info1]]),size = 2)+
              theme_light()+
              scale_color_gradient(low = "blue", high = "red") +
              xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))
          }else{ # if the color annotation is not a continuous variable
            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(colour = input[[group_info1]]),size = 2)+
              theme_light()+
              xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))+
              labs(colour="")
            
          }
          
          
        }
       else{ # if annotate by both color and shape
          if(is.numeric(scores[[input[[group_info1]]]])){

            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(shape = input[[group_info2]],colour = input[[group_info1]]),size = 2)+
              theme_light()+
              scale_color_gradient(low = "blue", high = "red") +
              xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))+
              scale_shape_manual(values=1:200)+
              labs(shape="")
          }else{

            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(shape = input[[group_info2]],colour = input[[group_info1]]),size = 2)+
              theme_light()+
              xlab(paste0("PC1"," (VarExp:",round(nipals_result$R2[1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(nipals_result$R2[2],4)*100,"%)"))+
              scale_shape_manual(values=1:200)+
              labs(colour="",shape="")
            
          }
          p

        }
        
        
      }

      
 #   }
    
 }
 , error=function(e){
   p=ggplot() + theme_void()+ggtitle("A continuous variable cannot be mapped to shape/It needs to have 2 or more levels for Voom.")
   p
 }
 )
}
