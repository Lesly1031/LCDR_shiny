# Trick file date creation update
onStop(function() {
  
  # File name
  p <- paste0("/srv/shiny-server/apps/production", "/app.R")
  
  # Update file 'date creation'
  Sys.setFileTime(p, now())
  
}) # onStop

study_id<-reactive({
  query = parseQueryString(session$clientData$url_search)
  #query$studyid
  if (!is.null(query[['cicpt_1']])) {
    cicpt_1 = query[["cicpt_1"]]
  } else{
    cicpt_1 = 1
  } 
  if (!is.null(query[['cicpt_2']])) {
    cicpt_2 = query[["cicpt_2"]]
  } else{
    cicpt_2=NULL
  }
  
  if(is.null(cicpt_2)){
    c(cicpt_1)
  } else{
    c(cicpt_1,cicpt_2)
  }
})


study_type<-reactive({
  query = parseQueryString(session$clientData$url_search)
  #query$studyid
  if (!is.null(query[['platform_1']])) {
    platform_1 = query[["platform_1"]]
  } else{
    platform_1 = 1
  } 
  if (!is.null(query[['platform_2']])) {
    platform_2 = query[["platform_2"]]
  } else{
    platform_2=NULL
  }
  
  if(is.null(platform_2)){
    c(platform_1)
  } else{
    c(platform_1,platform_2)
  }
})
########################### Boxplot UI########################
boxplot_choice = function(phenotype){
  phenotype = data.frame(phenotype)
  which_conti = sapply(names(phenotype),function(x) !is.numeric(phenotype[[x]]) )
  phenotype = data.frame(phenotype[names(which(which_conti))])
  names(phenotype)[-1]
}

box_plotUI<-function(headings,output_info,input_info,group_info,download_info,phenotype){
  
  column(8,
         
         panel(
           status="info",heading=paste0("Study: ",headings),
           dropdown(
             # autocomplete_input(input_info, label="Choose a gene",
             #                    options =geneID,value =geneID[1],max_options=20 ),
             
             selectizeInput(input_info, label="Choose a gene",
                            choices =NULL),
             # selectizeInput(input_info, label="Choose a gene",
             #                choices=NULL,selected =NULL),
             selectizeInput(group_info, label="Select group annotation",
                            choices=c(boxplot_choice(phenotype))),
             # selectizeInput(group_info, label="Select group annotation", 
             #                choices=NULL,selected =NULL),
             
             circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
             tooltip = tooltipOptions(title = "Click to see inputs !")
           ),
           conditionalPanel( 
             condition="input.update_box==0",
             br(),
             p("Please click 'Run' button to generate figures")
           ),
           
           conditionalPanel( 
             condition="input.update_box>0",
             br(),
             withLoader(plotlyOutput(output_info),type="html",loader="loader5")
             
           )
           
           
         )
  )
  
  
}

boxplot_observe = function(session,input_info,geneID,geneID_lip,study_type){
  if(is.null(geneID) || is.null(geneID_lip)){
    return()  # Exit if data is not available
  }
  if(grepl("lipidomics",study_type)){
    updateSelectizeInput(session=session, input_info, choices = list(
      "Collapse by class" = geneID_lip, 
      "Individual lipidomics" = geneID
    ))
    
  }else{
    updateSelectizeInput(session,input_info,choices = geneID,server = T)
    
  }
  
  
  
}


#-------------------------------------------------------------------------------
############### boxplot server with removing outlier black dots ##############
#-------------------------------------------------------------------------------
box_plotServer <- function(input,output,session,input_gene,input_group,data,group_data,folder_direct,study_type,data_lip,data_rnaseq){
  tryCatch({
    ########## change y labs of RNAseq, metabolomics and proteomics #############
    if(grepl("metabolomics",folder_direct)){
      lab_name = "Normalized Intensity"
    }else if (grepl("proteomics",folder_direct)){
      lab_name = "log2 intensity"
    }else{
      lab_name = "Voom normalized"
    }
    
    if(grepl("metabolomics",folder_direct)){
      lab_name_all = "Normalized Intensity"
    }else if (grepl("proteomics",folder_direct)){
      lab_name_all = "log2 intensity"
    }else{
      lab_name_all = "log2 TPM"
    }
    
    if(grepl("lipidomics",study_type)){
      common_cols <- intersect(names(data), names(data_lip))
      data <- data.frame(rbind(data[, common_cols], data_lip[, common_cols]),check.names=F)
      
    }else if(grepl("metabolomics",folder_direct)|grepl("proteomics",folder_direct)|grepl("lipidomics",folder_direct)){ 
      data=data.frame(data,check.names = F)
      
    }else{
      data = data.frame(data_rnaseq,check.names = F)
    }
    

    if(input[[input_group]]=="All"){
      
      data=data.frame(data,check.names = F)
      
      #}
      names_data=colnames(data)
      data=data.frame(data,check.names = F)
      names(data)=names_data
      group=data.frame(group_data)
      row_names_save=data$ID
      data=data.frame(t(data[,-c(1)]),check.names = F)
      names(data)=row_names_save
      
      data = data.frame(data[c(which(rownames(data)%in% group_data$X)),],check.names = F)
      
      
      for(i in 2:ncol(group)){
        group[,i]=as.character(group[,i])
      }
      #dt=data.frame(cbind(rownames(data),data[,"neg_00005_L.Alanine.D.Alanine.Beta.Alanine"]))
      dt=data.frame(cbind(rownames(data),data[,input[[input_gene]]]),check.names = F)
      
      
      dt[,2] = as.numeric(dt[,2])
      
      #names(dt)=c("X","neg_00005_L.Alanine.D.Alanine.Beta.Alanine")
      names(dt)=c("X",input[[input_gene]])
      
      p = ggplot(dt,aes_string(y=dt[,input[[input_gene]]]))+
        geom_boxplot()+
        # geom_jitter(aes_string(x=input[[input_group]],y=dt[,input[[input_gene]]],color=input[[input_group]]),shape=16, position=position_jitter(0.2))+
        theme_classic()+
        labs(x=input[[input_group]],y=lab_name,color=input[[input_group]])
      # ggtitle(paste0(input[[input_gene]]))
      
      # Add ggtitle based on input names
      gtitle = str_split(input[[input_gene]], "_")[[1]][1]
      if(gtitle=="X"){
        title_gg = input[[input_gene]]
      }else{
        title_gg = gtitle
      }
      
      p = p + ggtitle(title_gg)
      
      
      
      
      
      # p = ggplot(dt)+
      #   geom_boxplot(aes_string(y="neg_00005_L.Alanine.D.Alanine.Beta.Alanine"))+
      #   theme_classic()
      # 
      p
    }else {
      if(grepl("metabolomics",folder_direct)|grepl("proteomics",folder_direct)|grepl("lipidomics",folder_direct)){ ########### metabolomics ##################
        data=data.frame(data,check.names = F)
        
        #}
        names_data=colnames(data)
        data=data.frame(data,check.names = F)
        names(data)=names_data
        group=data.frame(group_data)
        row_names_save=data$ID
        data=data.frame(t(data[,-c(1)]),check.names = F)
        names(data)=row_names_save
        
        data = data.frame(data[c(which(rownames(data)%in% group$X)),],check.names = F)
        
        
        for(i in 2:ncol(group)){
          group[,i]=as.character(group[,i])
        }
        
        #dt=data.frame(cbind(rownames(data),data[,"ENSMUSG00000118380.1_AC131675.3"]))
        dt=data.frame(cbind(rownames(data),data[,input[[input_gene]]]))
        #clin = group[,c("X","histology")]
        clin = group[,c("X",input[[input_group]])]
        #clin$X = make.names(clin$X)
        dt = merge(dt,clin,by.x="X1",by.y="X")
        dt[,2] = as.numeric(dt[,2])
        
        #names(dt)=c("X","neg_00005_L.Alanine.D.Alanine.Beta.Alanine","histology")
        names(dt)=c("X",input[[input_gene]],input[[input_group]])
        
        p = ggplot(dt,aes_string(x=input[[input_group]],y=dt[[input[[input_gene]]]],color=input[[input_group]],label="X"))+
          geom_boxplot()+
          geom_jitter(shape=16, position=position_jitter(0.2))+
          theme_classic()+
          labs(x=input[[input_group]],y=lab_name,color=input[[input_group]])
        
        gtitle = str_split(input[[input_gene]], "_")[[1]][1]
        if(gtitle=="X"){
          title_gg = input[[input_gene]]
        }else{
          title_gg = gtitle
        }
        
        p = p + ggtitle(title_gg)
        
        p
        # p = ggplot(dt,aes_string(label= "X" ))+
        #   geom_boxplot(aes_string(x="histology",y="neg_00005_L.Alanine.D.Alanine.Beta.Alanine",color="histology"))+
        #   geom_jitter(aes_string(x="histology",y="neg_00005_L.Alanine.D.Alanine.Beta.Alanine",color="histology"),shape=16, position=position_jitter(0.2))+
        #   theme_classic()
        # 
        # p
      } else{ # RNAseq
        #### generate voom matrix ###
        dt = data.frame(data,check.names = F)
        dt_class=data.frame(group_data)
        
        ID_new = sapply(1:length(dt$ID),function(x){
          str_split(dt$ID[x], "_")[[1]][1]
        })
        dt = dt[which(!is.na(ID_new)),]
        
        dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names = F)
        
        rownames_dt_class=as.character(dt_class[,1])
        colnames_dt_class=colnames(dt_class)[-1]
        dt_class=data.frame(dt_class[,-1])
        rownames(dt_class)=rownames_dt_class
        names(dt_class)=colnames_dt_class
        
        group_info1 = input[[input_group]]
        
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
        
        dt_class = dt_class[match(names(dt_sub),rownames(dt_class)),]
        
        x <- DGEList(dt_sub) 
        grp  <- factor(dt_class[,group_info1])
        design <- model.matrix(~ 0+grp)
        colnames(design) = make.names(colnames(design))
        
        #keep.exprs <- filterByExpr(x, group = design)
        keep.exprs<-apply(dt_sub,1,sum)>=10 & rowSums(dt_sub!=0)>=4
        
        x <- x[keep.exprs, , keep.lib.sizes = FALSE]
        x <- calcNormFactors(x,method="TMM") # calculate normalization factor
        rnaseq_voomed = voom(x, design) 
        rnaseq_voomed = data.frame(t(rnaseq_voomed$E),check.names = F)
        
        rnaseq_voomed = merge(dt_class,rnaseq_voomed,by="row.names")
        
        
        p = ggplot(rnaseq_voomed,aes_string(x=input[[input_group]],y=input[[input_gene]],color=input[[input_group]],label="Row.names"))+
          geom_boxplot()+
          geom_jitter(shape=16, position=position_jitter(0.2))+
          theme_classic()+
          labs(x=input[[input_group]],y=lab_name,color=input[[input_group]])
        
        gtitle = str_split(input[[input_gene]], "_")[[1]][1]
        if(gtitle=="X"){
          title_gg = input[[input_gene]]
        }else{
          title_gg = gtitle
        }
        
        p = p + ggtitle(title_gg)
        
        p
      }
      
      
    }
    #p
    
    figure = ggplotly(p)%>% 
      layout(legend = list(orientation = 'h', x = 1.02, y = 1))
    
    for(i in 1:length(figure$x$data)){
      try({
        #only change geom_boxplot / box traces
        if (figure$x$data[[i]]$type == "box") {
          #set outliers to fully transparent
          figure$x$data[[i]]$marker$opacity = 0
        }
      })
    }
    
    figure
    
    
    
  },error = function(e){
    p= ggplot() + theme_void()+
      ggtitle("Only one group in the selected variable")
    ggplotly(p)
  })
  
  
}

#############################################################################
# voom download from deg
############################################################################
voom = fread("/voom_Study3990_rna-seq_flores_3990_rnaseq (1).csv")
voom
voom[1:5,1:5]
voom = data.frame(voom)
rownames(voom) = voom$V1;voom = data.frame(voom[,-1])
voom = data.frame(t(voom))

pheno =  fread("/RNAseq/phenotype_data.csv")
pheno = data.frame(pheno)
dt = merge(pheno,voom,by.x="X",by.y="row.names")

dt$EGFP_EGFP

ggplot(dt,aes_string(x="experimental_group",y=dt[["EGFP_EGFP"]],color="experimental_group",label="X"))+
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_classic()
