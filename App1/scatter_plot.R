##############Scatterplot function###############
##############Scatterplot function###############
scatter_plotUI<-function(headings,output_info,input_info1,input_info2,download_info,group_info,phenotype,spline_type){
  column(6,
         panel(
           status="info",heading=paste0("Study: ",headings),
           dropdown(
             # autocomplete_input(input_info1, label="Choose a gene", 
             #                    options =geneID,value =geneID[1],max_options=20),
             # autocomplete_input(input_info2, label="Choose another gene", 
             #                    options=geneID,value =geneID[2],max_options=20 ),
             # selectizeInput(input_info1, label="Choose a gene", 
             #                choices =geneID,selected =geneID[1]),
             # selectizeInput(input_info2, label="Choose another gene", 
             #                choices=geneID,selected =geneID[2]),
             
             selectizeInput(input_info1, label="Choose a gene", 
                            choices =NULL,selected =NULL),
             selectizeInput(input_info2, label="Choose another gene", 
                            choices=NULL,selected =NULL),
             selectizeInput(group_info, label="Select group annotation", 
                            choices =c("All",colnames(phenotype)[-1])),
             
            # radioButtons(spline_type,label="Add regression lines",choices=c("non-linear","linear"),selected = "linear"),
             downloadButton(download_info, "Download"),
             circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
             tooltip = tooltipOptions(title = "Click to see inputs !")
           ),
           
           conditionalPanel( 
             condition="input.update_scatter==0",
             br(),
             p("Please click 'Run' button to generate figures")
           ),
           
           conditionalPanel( 
             condition="input.update_scatter>0",
             br(),
             plotlyOutput(output_info)
           )
           

         )
  )
}

scatter_observe = function(session,input_info1,input_info2,geneID,geneID_lip,study_type){
    if(is.null(geneID) || is.null(geneID_lip)){
      return()  # Exit if data is not available
    }
  if(grepl("lipidomics",study_type)){
    updateSelectizeInput(session=session, input_info1, choices = list(
      "Collapse by class" = geneID_lip, 
      "Individual lipidomics" = geneID
    ),server=T)
    updateSelectizeInput(session=session, input_info2, choices = list(
      "Collapse by class" = geneID_lip, 
      "Individual lipidomics" = geneID
    ),selected =geneID_lip[2],server=T)
    
  }else{
    updateSelectizeInput(session,input_info1,choices = geneID,server=T)
    updateSelectizeInput(session,input_info2,choices = geneID,selected =geneID[2],server=T)
  }


}

scatter_plot <- function(input,output,session,data,group_data,input_group,input_gene1,input_gene2,study_type,data_lip){
  # if(input$scatter_meta=="Yes"){
  #   dt=data.frame(data_identify)
  #   
  # }else{
  #   dt=data.frame(data)
  #   
  # }
  if(grepl("lipidomics",study_type)){
    common_cols <- intersect(names(data), names(data_lip))
    data <- data.frame(rbind(data[, common_cols], data_lip[, common_cols]),check.names=F)
    
  }else{
    data=data.frame(data,check.names = F)
    
  }
  

  group_data=data.frame(group_data,check.names=F)
  data = data.frame(data[,c(1,which(names(data)%in% group_data$X))],check.names=F)
  
  names_data=colnames(data)
  data=data.frame(data,check.names=F)
  # for(i in 2:ncol(data)){
  #   data[,i] = log2(data[,i]+1)
  # }
  
  names(data)=names_data
  group=data.frame(group_data,check.names=F)
  row_names_save=data$ID
  data=data.frame(t(data[,-c(1)]),check.names=F)
  names(data)=row_names_save
  

  
  # for(i in 2:ncol(group)){
  #   group[,i]=as.character(group[,i])
  # }
  
  
  if(input[[input_group]]=="All"){
    
    
    dt=data.frame(cbind(rownames(data),data[,input[[input_gene1]]],data[,input[[input_gene2]]]),check.names=F)
    names(dt)=c("X",input[[input_gene1]],input[[input_gene2]])
    for(i in 2:3){
      dt[,i] = as.numeric(dt[,i])
    }
    t_pearson = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="pearson")$estimate),3)
    t_sp = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="spearman")$estimate),3)
    
    p = ggplot(dt,aes_string(test="X"))+
      geom_point(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]]))+
      theme_classic()+
      xlab(input[[input_gene1]])+
      ylab(input[[input_gene2]])+
      ggtitle(paste0("Pearson cor:",t_pearson,";Spearman cor:",t_sp))
    # if(input[[spline_type]]=="linear"){
    #   p= p+ geom_smooth(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]]),se=FALSE,method="lm")
    # }else if(input[[spline_type]]=="non-linear"){
    #   p= p+
    #     geom_smooth(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]]),se=FALSE,method="loess") 
    # }else{
    #   p =p
    # }
    
    p
    
  } 
  # else if(input[[input_group]]=="eligible samples (for metabolomics)"){
  #   #dt=data.frame(cbind(rownames(data),data[,"neg_00005_L.Alanine.D.Alanine.Beta.Alanine"],data[,"neg_00004_Glutathione"]))
  #   dt=data.frame(cbind(rownames(data),data[,input[[input_gene1]]],data[,input[[input_gene2]]]))
  #   #clin = group[,c("X","experimental_group")]
  #   clin = group
  #   clin$X = make.names(clin$X)
  #   dt = merge(dt,clin,by.x="X1",by.y="X")
  #   for(i in 2:3){
  #     dt[,i] = as.numeric(dt[,i])
  #   }
  #   #names(dt)=c("X","neg_00005_L.Alanine.D.Alanine.Beta.Alanine","neg_00004_Glutathione","experimental_group")
  #   names(dt)[1:3]=c("X","X1","X2")
  #   
  #   
  # 
  #   t_pearson = round(as.numeric(cor.test(dt$X1,dt$X2,method="pearson")$estimate),3)
  #   t_sp = round(as.numeric(cor.test(dt$X1,dt$X2,method="spearman")$estimate),3)
  #   
  #   p = ggplot(dt,aes_string(test="X"))+
  #     geom_point(aes(X1,X2))+
  #     theme_classic()+
  #     xlab(input[[input_gene1]])+
  #     ylab(input[[input_gene2]])+
  #     ggtitle(paste0("Pearson cor:",t_pearson,";Spearman cor:",t_sp))
  #   
  #   # if(input[[spline_type]]=="linear"){
  #   #   p= p+ geom_smooth(aes(X1,X2),se=FALSE,method="lm")
  #   # }else if(input[[spline_type]]=="non-linear"){
  #   #   p= p+
  #   #     geom_smooth(aes(X1,X2),se=FALSE,method="loess") 
  #   # }else{
  #   #   p = p
  #   # }
  # }
  else {
    
    #dt=data.frame(cbind(rownames(data),data[,"neg_00005_L.Alanine.D.Alanine.Beta.Alanine"],data[,"neg_00004_Glutathione"]))
    dt=data.frame(cbind(rownames(data),data[,input[[input_gene1]]],data[,input[[input_gene2]]]))
    #clin = group[,c("X","experimental_group")]
    clin = group[,c("X",input[[input_group]])]
   # clin$X = make.names(clin$X)
    dt = merge(dt,clin,by.x="X1",by.y="X")
    for(i in 2:3){
      dt[,i] = as.numeric(dt[,i])
    }
    #names(dt)=c("X","neg_00005_L.Alanine.D.Alanine.Beta.Alanine","neg_00004_Glutathione","experimental_group")
    names(dt)=c("X",input[[input_gene1]],input[[input_gene2]],input[[input_group]])
    
    
    t_pearson = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="pearson")$estimate),3)
    t_sp = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="spearman")$estimate),3)
    
    p = ggplot(dt,aes_string(test="X"))+
      geom_point(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]],color=input[[input_group]]))+
      theme_classic()+
      xlab(input[[input_gene1]])+
      ylab(input[[input_gene2]])+
      ggtitle(paste0("Pearson cor:",t_pearson,";Spearman cor:",t_sp))
    
    # if(input[[spline_type]]=="linear"){
    #   p= p+ geom_smooth(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]]),se=FALSE,method="lm")
    # }else if(input[[spline_type]]=="non-linear"){
    #   p= p+
    #     geom_smooth(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]]),se=FALSE,method="loess") 
    # }else{
    #   p=p
    # }
    
    p
    
    
  }  
  
}

