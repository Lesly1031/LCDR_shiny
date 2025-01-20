########################### Boxplot UI########################
boxplot_choice = function(phenotype){
  phenotype = data.frame(phenotype)
  which_conti = sapply(names(phenotype),function(x) !is.numeric(phenotype[[x]]) )
  phenotype = data.frame(phenotype[names(which(which_conti))])
  names(phenotype)[-1]
}


box_plotUI<-function(headings,output_info,input_info,group_info,phenotype){
  
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

boxplot_observe = function(session,input_info,geneID,study_type){
  if(is.null(geneID)){
    return()  # Exit if data is not available
  }

    updateSelectizeInput(session,input_info,choices = geneID,server = T)
  
}


#-------------------------------------------------------------------------------
############### boxplot server with removing outlier black dots ##############
#-------------------------------------------------------------------------------
box_plotServer <- function(input,output,session,input_gene,input_group,data,group_data){
  tryCatch({
    ########## change y labs of RNAseq, metabolomics and proteomics #############
      lab_name = "Normalized Intensity"

      lab_name_all = "Normalized Intensity"
      data=data.frame(data,check.names = F)
      
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


