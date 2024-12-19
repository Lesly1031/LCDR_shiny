immuneDeconvUI_heat<-function(headings,output_info,input_info,group_info,phenotype){
  column(10,
         panel(
           status="info",heading=paste0("Study: ",headings),
           dropdown(
             selectizeInput(group_info, label="Select group", 
                            choices =c("All",heat_choice(phenotype))),

             circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
             tooltip = tooltipOptions(title = "Click to see inputs !")
           ),
           
           conditionalPanel( 
             condition="input.update_heat_immuneDeconv==0",
             br(),
             p("Please generate the immune Deconvolution table prior to generating the heatmap. To generate the figures, please click on the 'Run' button.")
           ),
           
           conditionalPanel( 
             condition="input.update_heat_immuneDeconv>0",
             br(),
             (div(style='max-width: 100%; height: 500px; width: 90%;overflow-y: scroll;overflow-x: scroll;',
                   div(style="text-align: center;",
                      withLoader(plotOutput(output_info,width="100%",height="500px"),type="html",loader="loader5"))


             ))
            # DTOutput(output_info)
           )
           
           
           
         )
  )
  
}



####################### Use pheatmap ################################################
immuneDeconv_heat_server<-function(input,output,session,data,data_class,input_gene,group_info){
  
  
  tryCatch({
    
    # if(input$heat_meta=="Yes"){
    #   dt=data.frame(data_identify)
    #   
    # }else{
    #   dt=data.frame(data)
    #   
    # }
    #res = immunedeconv::deconvolute(data, "xcell")
    res = data # result of immune deconvolution
    dt=data.frame(res,check.names = F)
    dt_class = data.frame(data_class)
    
    
    dt_class$X[dt_class$X==""]=NA
    dt_class = dt_class[complete.cases(dt_class$X),]
    dt = data.frame(dt[,c(1,which(names(dt)%in% dt_class$X))],check.names = F)
    

    dt_sub = data.frame(dt)
    
    cell_type = dt_sub$cell_type
    

    dt_sub=data.frame(t(dt_sub[,-c(1)]),check.names=F)
    colnames(dt_sub) = cell_type


    cal_z_score <- function(x){
      (x - mean(x,na.rm=T)) / sd(x,na.rm = T)
    }

    dt_sub=apply(dt_sub, 2, cal_z_score)



    if(input[[group_info]]=="All"){
      dt_sub = data.frame(t(dt_sub),check.names = F)

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
      # length_nchar2 =  sapply(1:ncol(dt_sub),function(x){
      #   nchar(colnames(dt_sub)[x])
      # })
      #
      # length_nchar_l2 = which(length_nchar2>20)
      # colnames(dt_sub)[length_nchar_l2] = make.names(paste0(substr(colnames(dt_sub)[length_nchar_l2],1,20),"..."),unique=TRUE)

      if(nrow(dt_sub)==1){

        pheatmap(dt_sub[,1:(ncol(dt_sub))],cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                 show_rownames =T, show_colnames = T,main="z-score heatmap",na_col = "gray")
      }else{
        p= pheatmap(dt_sub[,1:(ncol(dt_sub))],cluster_rows=T,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                    show_rownames =T, show_colnames = T,main="z-score heatmap",na_col = "gray")

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

      dt_sub=data.frame(t(dt_sub),check.names = F)

      font_size=ifelse(nrow(dt_sub)<=20,12,ifelse(nrow(dt_sub)>20&nrow(dt_sub)<=50,10,
                                                  ifelse(nrow(dt_sub)>50&nrow(dt_sub)<=80,9,
                                                         ifelse(nrow(dt_sub)>80&nrow(dt_sub)<=100,8,
                                                                ifelse(nrow(dt_sub)>10&nrow(dt_sub)<=130,7,5)))))


      if(nrow(dt_sub)==1){

        pheatmap(dt_sub[,1:(ncol(dt_sub))],annotation = annotation,cluster_rows=F,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                 show_rownames =T, show_colnames = T,main="z-score heatmap",na_col = "gray")
      }else{
        p =  pheatmap(dt_sub[,1:(ncol(dt_sub))],annotation = annotation,cluster_rows=T,col = colorRampPalette(c("navy", "white", "firebrick3"))(50),fontsize=font_size,
                      show_rownames =T, show_colnames = T,main="z-score heatmap",na_col = "gray")
        p

      }
    }
    
  }

  )
}


