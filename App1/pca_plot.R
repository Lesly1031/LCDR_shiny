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


pca_UI<-function(headings,output_info,group_info1,group_info2,download_info,phenotype,identify_meta,study_type){
  if(grepl("metabolomics",study_type)){
    column(10,
           panel(
             status="info",heading=paste0("Study: ",headings),
             dropdownButton(
               # selectizeInput(exclude_info, label="Ignore unidentified metabolomics? (Not available for RNASeq)", 
               #                choices=c("Yes","No")),
               
               radioButtons(identify_meta,"Only include identified metabolites?",c("Yes","No"),"Yes"),
               
               selectizeInput(group_info1, label="Select color annotation", 
                              choices=c(names(phenotype)[-1]),selected = "experimental_group"),   
               selectizeInput(group_info2, label="Select shape annotation", 
                              choices=c("None",pca_choice_shape(phenotype)),selected = "None"),
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
  }else{
    column(10,
           panel(
             status="info",heading=paste0("Study: ",headings),
             dropdownButton(
               radioButtons(identify_meta,"",c(paste0("This is ",study_type," study")),paste0("This is ",study_type," study")),
               selectizeInput(group_info1, label="Select color annotation", 
                              choices=c(names(phenotype)[-1]),selected = "experimental_group"),   
               selectizeInput(group_info2, label="Select shape annotation", 
                              choices=c("None",pca_choice_shape(phenotype)),selected = "None"),    
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
  
  
  
}

#####################################PCA functions###############################################################
##############################PCA functions###################################################

pca_server<-function(input,output,session,data_identify,data,data_class,data_rna,group_info1,group_info2,folder_direct,identify_meta){
  
  if(input[[identify_meta]]=="Yes"){ 
    dt=data.frame(data_identify,check.names = F)
    
  }else{
    if(grepl("metabolomics",folder_direct)|grepl("proteomics",folder_direct)|grepl("lipidomics",folder_direct)){
      dt=data.frame(data,check.names = F)
      
    }else{
      dt=data.frame(data_rna,check.names = F)
    }

  }
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

      if(grepl("metabolomics",folder_direct)|grepl("proteomics",folder_direct)|grepl("lipidomics",folder_direct)){
        dt_sub = dt_sub[ , which(apply(dt_sub, 2, var,na.rm=T) != 0)]
        
        nipals_result <- nipals(dt_sub, scale= TRUE, center=TRUE,ncomp = 2)
        scores = data.frame(nipals_result$scores,check.names = F)
        scores$X = rownames(scores)
        scores = merge(scores,dt_class,by.x="X",by.y="row.names")
        for(i in 2:3){
          scores[,i]=as.numeric(as.character(scores[,i]))
        }
        
 
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
      else{
        
        names_na = rownames(dt_class)[which(is.na(dt_class[input[[group_info1]]]))]
        dt_sub = data.frame(t(dt_sub),check.names=F)
        
        
        if (length(names_na) != 0) {
          dt_sub = dt_sub[, -(which(colnames(dt_sub) %in% names_na))]
          dt_class = na.omit(dt_class)
        } else {
          dt_sub = dt_sub
          dt_class = dt_class
        }
        dt_class = dt_class[match(names(dt_sub), rownames(dt_class)), ]
        x <- DGEList(dt_sub)
        
        low_expressed_filtering <- function(x, N_min) {
          # Genes with 0 count for all samples will be removed N_min == minimum
          # number of samples without 0 count; it won't change PCA plot much, will
          # change the number of DEG list
          keep.exprs <- apply(dt_sub, 1, sum) > 0 & rowSums(dt_sub != 0) >= N_min
          x <<- x[keep.exprs, , keep.lib.sizes = FALSE]
          ## voom(x,plot=T)
        }
        
        low_expressed_filtering(x, 20)
        x <- calcNormFactors(x, method = "TMM")  # calculate normalization factor: according to Chia-Ho's email on 01/11/2023
        
        grp <- factor(dt_class[[input[[group_info1]]]])
        design <- model.matrix(~0 + grp)
        colnames(design) = make.names(colnames(design))
        
        rnaseq_voomed = voom(x, design)
        rnaseq_voomed_out <- as.data.frame(rnaseq_voomed$E)
        
        pca_res <- prcomp(t(rnaseq_voomed_out), center = T, scale. = T)
        t = data.frame(pca_res$x, check.names = F)
        scores = data.frame(X = rownames(t), PC1 = t$PC1, PC2 = t$PC2, check.names = F)
        scores = merge(scores, dt_class, by.x = "X", by.y = "row.names")
        
        for (i in 2:3) {
          scores[, i] = as.numeric(as.character(scores[, i]))
        }


       if(input[[group_info2]]=="None"){ # if annotate by color
          if(is.numeric(scores[[input[[group_info1]]]])){ # if the color annotation is a continuous variable
            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(colour = input[[group_info1]]),size = 2)+
              theme_light()+
              scale_color_gradient(low = "blue", high = "red") +
              xlab(paste0("PC1"," (VarExp:",round(summary(pca_res)$importance[2,1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(summary(pca_res)$importance[2,2],4)*100,"%)"))
          }else{ # if the color annotation is not a continuous variable
            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(colour = input[[group_info1]]),size = 2)+
              theme_light()+
              xlab(paste0("PC1"," (VarExp:",round(summary(pca_res)$importance[2,1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(summary(pca_res)$importance[2,2],4)*100,"%)"))+
              labs(color = input[[group_info1]])
          }
          p
          
        }
        else{ # if annotate by both color and shape
          if(is.numeric(scores[[input[[group_info1]]]])){
            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(shape = input[[group_info2]],colour = input[[group_info1]]),size = 2)+
              theme_light()+
              scale_color_gradient(low = "blue", high = "red") +
              xlab(paste0("PC1"," (VarExp:",round(summary(pca_res)$importance[2,1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(summary(pca_res)$importance[2,2],4)*100,"%)"))+
              scale_shape_manual(values=1:200)+
              labs(shape="")
          }else{
            p= ggplot(scores,aes(x=PC1,y=PC2,text = X
            ))+
              geom_point(aes_string(colour = input[[group_info1]],shape = input[[group_info2]]),size = 2)+
              theme_light()+
              xlab(paste0("PC1"," (VarExp:",round(summary(pca_res)$importance[2,1],4)*100,"%)"))+
              ylab(paste0("PC2"," (VarExp:",round(summary(pca_res)$importance[2,2],4)*100,"%)"))+
              labs(color = input[[group_info1]],shape=input[[group_info2]])+
              scale_shape_manual(values=1:200)+
              ggtitle("Voom is based on the selection of color annotation")
          }
          
        }
        
        
        p

      }
      
 #   }
    
 }
 , error=function(e){
   p=ggplot() + theme_void()+ggtitle("A continuous variable cannot be mapped to shape/It needs to have 2 or more levels for Voom.")
   p
 }
 )
}
