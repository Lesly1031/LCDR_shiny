immuneDeconv_box_choice = function(phenotype){
  phenotype = data.frame(phenotype)
  which_conti = sapply(names(phenotype),function(x) !is.numeric(phenotype[[x]]) )
  phenotype = data.frame(phenotype[names(which(which_conti))])
  names(phenotype)[-1]
}

immuneDeconvUI_box<-function(headings,output_info,group_info,phenotype,cell_type){
  column(10,
         panel(
           status="info",heading=paste0("Study: ",headings),
           dropdown(
             selectInput(cell_type,"Choose a cell type", choices=NULL),
             selectizeInput(group_info, label="Select group", 
                            choices =c("All",immuneDeconv_box_choice(phenotype))),
             
             circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
             tooltip = tooltipOptions(title = "Click to see inputs !")
           ),
           
           conditionalPanel( 
             condition="input.update_box_immuneDeconv==0",
             br(),
             p("Please generate the immune Deconvolution table prior to generating the boxplot To generate the figures, please click on the 'Run' button.")
           ),
           
           conditionalPanel( 
             condition="input.update_box_immuneDeconv>0",
             br(),
             (div(style='max-width: 100%; height: 500px; width: 90%;overflow-y: scroll;overflow-x: scroll;',
                  div(style="text-align: center;",
                      withLoader(plotlyOutput(output_info,width="100%",height="500px"),type="html",loader="loader5"))
                  
                  
             ))
             # DTOutput(output_info)
           )
           
           
           
         )
  )
  
}

immuneDeconvUI_method_sel_box = function(session,data,cell_type){
  
  # data is results from immuneDeconv table
  
    updateSelectizeInput(session,cell_type,choices = data$cell_type ,server = T)

  
}

immuneDeconv_box_server <- function(input,output,session,input_gene,input_group,data,group_data,folder_direct){
# input_gene refers to input immune cell type
# data refers to the immunedeconv result
  tryCatch({
    if(input[[input_group]]=="All"){
      data=data.frame(data,check.names = F)
      
      names_data=colnames(data)
      data=data.frame(data,check.names = F)
      names(data)=names_data
      group=data.frame(group_data)
      row_names_save=data$cell_type
      data=data.frame(t(data[,-c(1)]),check.names = F)
      names(data)=row_names_save
      
      data = data.frame(data[c(which(rownames(data)%in% group_data$X)),],check.names = F)
      
      
      for(i in 2:ncol(group)){
        group[,i]=as.character(group[,i])
      }

      dt=data.frame(cbind(rownames(data),data[,input[[input_gene]]]),check.names = F)
      
      
      dt[,2] = as.numeric(dt[,2])
      
      names(dt)=c("X",input[[input_gene]])
      
      p = ggplot(dt,aes_string(y=dt[,input[[input_gene]]]))+
        geom_boxplot()+
        theme_classic()+
        labs(x=input[[input_group]],y="",color=input[[input_group]])

      # Add ggtitle based on input names
      p = p + ggtitle( input[[input_gene]])
    }else{
      data=data.frame(data,check.names = F)
      
      names_data=colnames(data)
      data=data.frame(data,check.names = F)
      names(data)=names_data
      group=data.frame(group_data)
      row_names_save=data$cell_type
      data=data.frame(t(data[,-c(1)]),check.names = F)
      names(data)=row_names_save
      
      data = data.frame(data[c(which(rownames(data)%in% group$X)),],check.names = F)
      
      
      for(i in 2:ncol(group)){
        group[,i]=as.character(group[,i])
      }
      
      dt=data.frame(cbind(rownames(data),data[,input[[input_gene]]]))
      clin = group[,c("X",input[[input_group]])]
      dt = merge(dt,clin,by.x="X1",by.y="X")
      dt[,2] = as.numeric(dt[,2])
      
      names(dt)=c("X",input[[input_gene]],input[[input_group]])
      
      p = ggplot(dt,aes_string(x=input[[input_group]],y=dt[[input[[input_gene]]]],color=input[[input_group]],label="X"))+
        geom_boxplot()+
        geom_jitter(shape=16, position=position_jitter(0.2))+
        theme_classic()+
        labs(x=input[[input_group]],y="",color=input[[input_group]])
      
      
      p = p + ggtitle(input[[input_gene]])
      
      p
    }
    
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
