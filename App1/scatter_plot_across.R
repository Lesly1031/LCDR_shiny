scatter_across_plotUI<-function(headings,output_info,input_info1,input_info2,download_info,group_info,phenotype){
  column(6,
         panel(
           status="info",heading=paste0("Study: ",headings),
           dropdown(
             
             selectizeInput(input_info1, label="Choose a gene from 1st study", 
                            choices =NULL,selected =NULL),
             selectizeInput(input_info2, label="Choose a gene from 2nd study", 
                            choices=NULL,selected =NULL),
             selectizeInput(group_info, label="Select group annotation", 
                            choices =c("All",colnames(phenotype)[-1])),
             
             downloadButton(download_info, "Download"),
             circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
             tooltip = tooltipOptions(title = "Click to see inputs !")
           ),
           
           conditionalPanel( 
             condition="input.update_scatter_across==0",
             br(),
             p("Please click 'Run' button to generate figures")
           ),
           
           conditionalPanel( 
             condition="input.update_scatter_across>0",
             br(),
             plotlyOutput(output_info)
             # DTOutput(output_info)
           )
           
         )
  )
}

scatter_observe_across = function(session,input_info1,input_info2,geneID1,geneID2){
  updateSelectizeInput(session,input_info1,choices = geneID1,server = T)
  updateSelectizeInput(session,input_info2,choices = geneID2,server = T)
  
}

scatter_plot_across <- function(input,output,session,data1,data2,group_data,input_group,input_gene1,input_gene2){
  
  data1 = data.frame(data1,check.names = F)
  data2 =  data.frame(data2,check.names = F) 
  pheno =  data.frame(group_data,check.names = F) 
  
  rownames(data1) = data1$ID;data1 = data.frame(data1[,-1],check.names = F)
  rownames(data2) = data2$ID;data2 = data.frame(data2[,-1],check.names = F)
  
  data1 = data.frame(t(data1),check.names = F)
  data2 = data.frame(t(data2),check.names = F)
  
  
  data1_sub = data.frame(data1[[input[[input_gene1]]]],row.names = rownames(data1),check.names=F)
  names(data1_sub) = input[[input_gene1]]
  data2_sub = data.frame(data2[[input[[input_gene2]]]],row.names = rownames(data2),check.names=F)
  names(data2_sub) = input[[input_gene2]]
  
  
  data = merge(data1_sub,data2_sub,by="row.names",check.names=F)
  
  group=data.frame(pheno,check.names=F)
  data = merge(data,group,by.x="Row.names",by.y="X",check.names=F)
  
  
  if(input[[input_group]]=="All"){
    
    
    dt=data.frame(cbind(data$Row.names,data[,input[[input_gene1]]],data[,input[[input_gene2]]]),check.names=F)
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
    
    p
    
  }  else {
    
    #dt=data.frame(cbind(rownames(data),data[,"neg_00005_L.Alanine.D.Alanine.Beta.Alanine"],data[,"neg_00004_Glutathione"]))
    dt=data.frame(cbind(data$Row.names,data[,input[[input_gene1]]],data[,input[[input_gene2]]]))
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
    
    
    
    p
    
    
  }
  
}

# scatter_across_plotUI<-function(headings,output_info,input_info1,input_info2,download_info,group_info,phenotype,lip_class1,lip_class2,study_type1,study_type2){
#   if(grepl("lipidomics",study_type1)){
#     column(6,
#            panel(
#              status="info",heading=paste0("Study: ",headings),
#              dropdown(
#                radioButtons(lip_class1,"Collapse class for lipidomics?",c("Yes","No"),"Yes"),
#                
#                selectizeInput(input_info1, label="Choose a gene from 1st study", 
#                               choices =NULL,selected =NULL),
#                selectizeInput(input_info2, label="Choose a gene from 2nd study", 
#                               choices=NULL,selected =NULL),
#                selectizeInput(group_info, label="Select group annotation", 
#                               choices =c("All",colnames(phenotype)[-1])),
#                
#                downloadButton(download_info, "Download"),
#                circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
#                tooltip = tooltipOptions(title = "Click to see inputs !")
#              ),
#              
#              conditionalPanel( 
#                condition="input.update_scatter_across==0",
#                br(),
#                p("Please click 'Run' button to generate figures")
#              ),
#              
#              conditionalPanel( 
#                condition="input.update_scatter_across>0",
#                br(),
#                plotlyOutput(output_info)
#                #  DTOutput(output_info)
#              )
#              
#            )
#     )
#   }else if(grepl("lipidomics",study_type2)){
#     column(6,
#            panel(
#              status="info",heading=paste0("Study: ",headings),
#              dropdown(
#                radioButtons(lip_class2,"Collapse class for lipidomics?",c("Yes","No"),"Yes"),
#                
#                selectizeInput(input_info1, label="Choose a gene from 1st study", 
#                               choices =NULL,selected =NULL),
#                selectizeInput(input_info2, label="Choose a gene from 2nd study", 
#                               choices=NULL,selected =NULL),
#                selectizeInput(group_info, label="Select group annotation", 
#                               choices =c("All",colnames(phenotype)[-1])),
#                
#                downloadButton(download_info, "Download"),
#                circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
#                tooltip = tooltipOptions(title = "Click to see inputs !")
#              ),
#              
#              conditionalPanel( 
#                condition="input.update_scatter_across==0",
#                br(),
#                p("Please click 'Run' button to generate figures")
#              ),
#              
#              conditionalPanel( 
#                condition="input.update_scatter_across>0",
#                br(),
#                plotlyOutput(output_info)
#                #  DTOutput(output_info)
#              )
#              
#            )
#     )
#   }else{
#     column(6,
#            panel(
#              status="info",heading=paste0("Study: ",headings),
#              dropdown(
#                
#                selectizeInput(input_info1, label="Choose a gene from 1st study", 
#                               choices =NULL,selected =NULL),
#                selectizeInput(input_info2, label="Choose a gene from 2nd study", 
#                               choices=NULL,selected =NULL),
#                selectizeInput(group_info, label="Select group annotation", 
#                               choices =c("All",colnames(phenotype)[-1])),
#                
#                downloadButton(download_info, "Download"),
#                circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
#                tooltip = tooltipOptions(title = "Click to see inputs !")
#              ),
#              
#              conditionalPanel( 
#                condition="input.update_scatter_across==0",
#                br(),
#                p("Please click 'Run' button to generate figures")
#              ),
#              
#              conditionalPanel( 
#                condition="input.update_scatter_across>0",
#                br(),
#                plotlyOutput(output_info)
#                #  DTOutput(output_info)
#              )
#              
#            )
#     )
#   }
# 
# }
# 
# scatter_observe_across = function(input,session,input_info1,input_info2,geneID1,geneID2,geneID1_class,geneID2_class,lip_class1,lip_class2){
#   
#   # if (!is.null(input[[lip_class]]) && length(input[[lip_class]]) > 0) {
#   #   if(input[[lip_class]] == "Yes") {
#   #     # Update with geneID_class if "Yes"
#   #     updateSelectizeInput(session,input_info1,choices = geneID1_class,server = T)
#   #     updateSelectizeInput(session,input_info2,choices = geneID2_class,server = T)
#   #   } else {
#   #     # Update with geneID if "No"
#   #     updateSelectizeInput(session,input_info1,choices = geneID1,server = T)
#   #     updateSelectizeInput(session,input_info2,choices = geneID2,server = T)
#   #   }
#   #   
#   # }
#   # Check for lip_class1 and update input_info1 accordingly
#   if (!is.null(input[[lip_class1]]) && length(input[[lip_class1]]) > 0) {
#     if(input[[lip_class1]] == "Yes") {
#       # Update with geneID1_class if "Yes"
#       updateSelectizeInput(session, input_info1, choices = geneID1_class, server = TRUE)
#     } else {
#       # Update with geneID1 if "No"
#       updateSelectizeInput(session, input_info1, choices = geneID1, server = TRUE)
#     }
#   }
#   
#   # Check for lip_class2 and update input_info2 accordingly
#   if (!is.null(input[[lip_class2]]) && length(input[[lip_class2]]) > 0) {
#     if(input[[lip_class2]] == "Yes") {
#       # Update with geneID2_class if "Yes"
#       updateSelectizeInput(session, input_info2, choices = geneID2_class, server = TRUE)
#     } else {
#       # Update with geneID2 if "No"
#       updateSelectizeInput(session, input_info2, choices = geneID2, server = TRUE)
#     }
#   }
#   
# 
#   
# }
# 
# scatter_plot_across <- function(input,output,session,data1,data2,group_data,input_group,input_gene1,input_gene2,study_type1,study_type2,lip_class1,lip_class2,data_lip1,data_lip2){
# 
#   if(input[[study_type1]]=="lipidomics"&&input[[lip_class1]]=="Yes"){
#     data1 = data_lip1
#   }else if(input[[study_type2]]=="lipidomics"&&input[[lip_class2]]=="Yes"){
#     data2 = data_lip2
#     
#   }
#   
#   data1 = data.frame(data1,check.names = F)
#   data2 =  data.frame(data2,check.names = F) 
#   pheno =  data.frame(group_data,check.names = F) 
#   
#   rownames(data1) = data1$ID;data1 = data.frame(data1[,-1],check.names = F)
#   rownames(data2) = data2$ID;data2 = data.frame(data2[,-1],check.names = F)
#   
#   data1 = data.frame(t(data1),check.names = F)
#   data2 = data.frame(t(data2),check.names = F)
#   
# 
#   data1_sub = data.frame(data1[[input[[input_gene1]]]],row.names = rownames(data1),check.names=F)
#   names(data1_sub) = input[[input_gene1]]
#   data2_sub = data.frame(data2[[input[[input_gene2]]]],row.names = rownames(data2),check.names=F)
#   names(data2_sub) = input[[input_gene2]]
#   
# 
#   data = merge(data1_sub,data2_sub,by="row.names",check.names=F)
# 
#   group=data.frame(pheno,check.names=F)
#   data = merge(data,group,by.x="Row.names",by.y="X",check.names=F)
#   
#   
#   if(input[[input_group]]=="All"){
# 
# 
#     dt=data.frame(cbind(data$Row.names,data[,input[[input_gene1]]],data[,input[[input_gene2]]]),check.names=F)
#     names(dt)=c("X",input[[input_gene1]],input[[input_gene2]])
#     for(i in 2:3){
#       dt[,i] = as.numeric(dt[,i])
#     }
#     t_pearson = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="pearson")$estimate),3)
#     t_sp = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="spearman")$estimate),3)
# 
#     p = ggplot(dt,aes_string(test="X"))+
#       geom_point(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]]))+
#       theme_classic()+
#       xlab(input[[input_gene1]])+
#       ylab(input[[input_gene2]])+
#       ggtitle(paste0("Pearson cor:",t_pearson,";Spearman cor:",t_sp))
# 
#     p
# 
#   }  else {
# 
#     #dt=data.frame(cbind(rownames(data),data[,"neg_00005_L.Alanine.D.Alanine.Beta.Alanine"],data[,"neg_00004_Glutathione"]))
#     dt=data.frame(cbind(data$Row.names,data[,input[[input_gene1]]],data[,input[[input_gene2]]]))
#     #clin = group[,c("X","experimental_group")]
#     clin = group[,c("X",input[[input_group]])]
#     # clin$X = make.names(clin$X)
#     dt = merge(dt,clin,by.x="X1",by.y="X")
#     for(i in 2:3){
#       dt[,i] = as.numeric(dt[,i])
#     }
#     #names(dt)=c("X","neg_00005_L.Alanine.D.Alanine.Beta.Alanine","neg_00004_Glutathione","experimental_group")
#     names(dt)=c("X",input[[input_gene1]],input[[input_gene2]],input[[input_group]])
# 
# 
#     t_pearson = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="pearson")$estimate),3)
#     t_sp = round(as.numeric(cor.test(dt[[input[[input_gene1]]]],dt[[input[[input_gene2]]]],method="spearman")$estimate),3)
# 
#     p = ggplot(dt,aes_string(test="X"))+
#       geom_point(aes_string(dt[,input[[input_gene1]]],dt[,input[[input_gene2]]],color=input[[input_group]]))+
#       theme_classic()+
#       xlab(input[[input_gene1]])+
#       ylab(input[[input_gene2]])+
#       ggtitle(paste0("Pearson cor:",t_pearson,";Spearman cor:",t_sp))
# 
# 
# 
#     p
# 
# 
#   }
#   
# }
# 
