library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ChIPseeker)
library(rtracklayer)
library(RColorBrewer)
library(VennDiagram)
library(fst)
library(shinycustomloader)
library(DT)
library(RColorBrewer)
library(data.table)
library(shiny)
library(shinydashboardPlus)
library(shinyWidgets)
library(limma)
library(edgeR)
library(stringr)
library(ggplot2)
library(plotly)
library(tidyverse)
library(caTools)
library(shinyjs)
library(reshape2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(rhandsontable)

source("/venn_diag.R")
source("/peak_distribution.R")
source("/peak_heat.R")
source("/peak_meta.R")
source("/peak_heat_tss.R")
source("/peak_tss_meta.R")
source("/DEG_table.R")
source("/DEG_fig.R")
source("/sample_table.R")

#load("/files/study/demo/region_plot/bwgr.RData")


# Define UI for application that draws a histogram
ui <- fluidPage(
  shinyjs::useShinyjs(),
  fluidRow(
    br(),
    tabsetPanel(id="mainnav",
      tabPanel("Phenotype",
               fluidRow(
                 column(2,
                        actionButton("instruction_pheno","Tab instructions")
                 ),
                 column(10,
                        withLoader( rHandsontableOutput("all_sample"),type="html",loader="loader5"),
                        tableOutput("test")
                 )
               )
               
      ),
 
      tabPanel("Peak distribution",
               fluidRow(
                 br(),
                 style = "margin-left: 20px;",
                 column(3,
                        selectizeInput("sample_sel_peak_distr","Select samples from the dropdown list",choices =NULL,multiple = T),
                        p("Click 'Run' button to get the figure"),
                        actionButton("update_peak_distribution","Run", icon("refresh"))
                 ),
                 column(9,
                        uiOutput("peak_distribution_ui")
                 )

               )),
      tabPanel("Peak heatmap",
               fluidRow(
                 style = "margin-left: 20px;",
                 column(2,
                        br(),
                        selectizeInput("sample_sel_peak_heat","Select samples from the dropdown list",choices =NULL,multiple = T,size = 2,options = list(maxItems = 4)),
                        p("Click 'Run' button to get the figure"),
                        actionButton("update_peak_heat","Run", icon("refresh"))
                 ),
                 column(10,
                      uiOutput("peak_heat")
                 )
                 
                 
               )
      ),
      tabPanel("Peak meta",
               fluidRow(
                 style = "margin-left: 20px;",
                 column(2,
                        br(),
                        selectizeInput("sample_sel_peak_meta","Select samples from the dropdown list",choices =NULL,multiple = T,size = 2,options = list(maxItems = 4)),
                        p("Click 'Run' button to get the figure"),
                        actionButton("update_peak_meta","Run", icon("refresh"))
                 ),
                 column(10,
                        uiOutput("peak_meta")
                 )
               )

      ),
      tabPanel("Peak TSS heatmap",
               fluidRow(
                 style = "margin-left: 20px;",
                 column(2,
                        br(),
                        selectizeInput("sample_sel_peak_TSS_heat","Select samples from the dropdown list",choices =NULL,multiple = T, size = 2,options = list(maxItems = 4)),
                        p("Click 'Run' button to get the figure"),
                        actionButton("update_peak_tss_heat","Run", icon("refresh"))
                 ),
                 column(10,
                        uiOutput("peak_tss_heat")
                        
                 )
               )
      ),
      tabPanel("Peak TSS meta",
               fluidRow(
                 style = "margin-left: 20px;",
                 column(2,
                        br(),
                        selectizeInput("sample_sel_peak_TSS_meta","Select samples from the dropdown list",choices =NULL,multiple = T, size = 2,options = list(maxItems = 4)),
                        p("Click 'Run' button to get the figure"),
                        actionButton("update_peak_tss_meta","Run", icon("refresh"))
                 ),
                 column(10,
                        uiOutput("peak_tss_meta")
                        
                 )
                 
                 
               )
      ),

      tabPanel("DEG",
               tabsetPanel(type="tabs",id="mainnav_deg",

                           tabPanel("DEG test result",
                                    fluidRow(
                                      style = "margin-left: 20px;",
                                      column(2,
                                             br(),
                                             selectizeInput("deg_group", label="Select group annotation", choices =NULL),
                                             selectizeInput("region_deg","Select Chr from the dropdown list",choices =c("All",paste0("chr",1:22),"chrX","chrY"),selected="All",
                                                            multiple = F),
                                             numericInput("start_deg","Start position",value=0,min=0),
                                             numericInput("end_deg","End position",value=0,min=0),
                                             actionButton("update_DEG","Run", icon("refresh")),
                                             downloadButton("download_DEG", "Download")
                                      ),
                                      column(10,
                                             conditionalPanel( 
                                               condition="input.update_DEG==0",
                                               br(),
                                               p("Please click 'Run' button to generate figures")
                                             ),
                                             
                                             conditionalPanel( 
                                               condition="input.update_DEG>0",
                                               br(),
                                               withLoader(DTOutput("deg_table"),type = "html",loader = "dnaspin")
                                               
                                             )
                                      )
                                      
                                    )       
                           ),
                           tabPanel("Vocalno plot and MA plot",
                                    fluidRow(
                                      style = "margin-left: 20px;",
                                      column(2,selectInput('plot_type_MA',"Select a figure type",choices = c("MA plot","Volcano plot")),
                                             selectInput("which_p","Choose p value or adjusted p",choices = c("pval","adjust.p"),selected="pval"),
                                             numericInput("p_threshold_MA", "pval/adj.pval threshold", 0.05,min=0.0001,max=0.1,step=0.05),
                                             numericInput("fc_threshold_MA", "Select fold change threshold", value=1.5),
                                             selectizeInput("deg_fill", 'Optionally, input interested genes to highlight',
                                                            choices =NULL,multiple=T),
                                             actionButton("update_MA","Run", icon("refresh"))
                                      ),
                                      column(10,
                                             conditionalPanel( 
                                               condition="input.update_MA==0",
                                               br(),
                                               p("Please click 'Run' button to generate figures")
                                             ),
                                             
                                             conditionalPanel( 
                                               condition="input.update_MA>0",
                                               br(),
                                              # DTOutput("deg_fig")
                                               withLoader(plotlyOutput("deg_fig",width = "50%",height = "300px"),type="html",loader="dnaspin")
                                             )
                                            )
                                    )  
                           ),
                           tabPanel(
                             "DEG Region plot",
                             fluidRow(
                               style = "margin-left: 20px;",
                               column(2,
                                      br(),
                                      selectizeInput("chr_region_deg","Select Chr from the dropdown list",choices =c(paste0("chr",1:22),"chrX","chrY"),
                                                     selected="chr4",multiple = F),
                                      numericInput("start_deg_region","Start position",value = 113853504, min=1),
                                      numericInput("end_deg_region","End position",value = 114476778, min=1),
                                      selectizeInput("gene_region_deg","Optionally, Select a gene from the dropdown list",choices=NULL),
                                      #numericInput("window_deg","Width of moving window",value = 20, min=1),
                                      selectizeInput("sel_sample_regionDeg","Select samples to show the region plot in one figure",choices=NULL,multiple=T,options = list(maxItems = 4)),
                                      p("Click 'Run' button to get the figure"),
                                      actionButton("update_region_deg","Run", icon("refresh"))
                               ),
                               column(10,
                                      fluidRow(
                                        column(10,
                                              withLoader(plotOutput("region_plot_deg"),type = "html",loader = "dnaspin"),
                                               )),
                                      fluidRow(
                                        column(10,
                                               plotOutput("region_plot_DEG_each",height = "500px")
                                        )
                                      ),
                                      fluidRow(
                                        column(10,
                                               plotOutput("region_plot_DEG_overlap")
                                        )
                                      )

                               )
                               
                               
                             )
                           )
                           
               )
               
      )
    
    )
    
  )
)


server <- function(input,output,session){
  #---------------------------------------------------
  ######## Files for preparation ###################
  #---------------------------------------------------
  file_folder = reactive({
    "/files/study/3675/"
  })
  file_venn = reactive({
    "/"
    
  })
  beds = reactive({
    beds = list.files(paste0(file_folder(),"beds/"),pattern = "\\.narrowPeak$")
    beds
  })
  
  bigwig = reactive({
    bigwig = list.files(paste0(file_folder(),"bigwig/"),pattern = "\\.bw$")
    bigwig
  })
  
  samples = reactive({
    sapply( strsplit( bigwig(), "\\."), "[", 1)
  })
  
  samples_bed = reactive({
    
    sapply( strsplit( beds(), "\\."), "[", 1)
  })
  
  

  #===============================================================================
  ########################## phenotype data file ##############################
  #===============================================================================
  onclick('instruction_pheno', showModal(modalDialog(size="l",fade=TRUE,
                                                     title = "",
                                                     renderUI(
                                                       tags$img(src="images/pheno_info.jpg",width="800px")
                                                     ),easyClose = TRUE
  )))
  

  # input update
  phenotype = reactive({
    fread(paste0(file_folder(),"phenotype.csv"))
  })
  
  
  
  sample_file = reactive({
    rhandsontable(phenotype(),width = 800, height = 500, stretchH = "all")
  })
  
  output$all_sample = renderRHandsontable({
    sample_file()
  })
  
  
  phenotype_to_use <- reactive({
    dt_class = hot_to_r(input[["all_sample"]])
    
    # Debug output right after conversion

    if(is.null(dt_class)) {
      return(data.frame(X=character()))
    }
    
    dt_class = data.frame(dt_class, check.names = FALSE)
    dt_class = dt_class[complete.cases(dt_class$X),]
    
    # Final output before it's used

    return(dt_class)
  })
  
  pheno_name = reactive({
    if(is.null(phenotype_to_use()) || ncol(phenotype_to_use()) == 0) {
      return(NULL)
    }
    res = phenotype_to_use()$X
    res[res == ""] = NA
    na.omit(res)
  })
  
  output$test = renderTable({
    res= phenotype_to_use()$X
    res[res==""] = NA
    na.omit(res)
  })

  
  
  #-----------------------------------------------------------------------------
  ######################## Venn diagram #####################################
  #-----------------------------------------------------------------------------
  # input update
  observe({
    if(is.null(samples_bed()) || is.null(pheno_name())) {
      return()
    }
    choices = samples_bed()[samples_bed() %in% pheno_name()]
    updateSelectizeInput(session, "sample_sel_venn", choices = choices, selected = choices[1:3])
  })
  
  
  # bed file import
  bedtab =  eventReactive(input$update_venn,{
    lapply(1:length(input$sample_sel_venn),function(x){
      fread(paste0(file_folder(),"/beds/",input$sample_sel_venn[x],".narrowPeak"))
    })
  }
  )
  bedgr = eventReactive(input$update_venn,{
    bedgr = sortSeqlevels(GRangesList(lapply(bedtab(),makeGRangesFromDataFrame,seqnames.field='V1',
                                             start.field='V2', end.field='V3',keep.extra.columns=T)))
    seqlengths(bedgr) = seqlengths(Mmusculus)[seqlevels(bedgr)]
    bedgr
  })

  venn_diagram = eventReactive(input$update_venn,{
    venn_gram(input,output,session,bedgr())

  })
  output$venn_diagram = renderPlot({
    venn_diagram()
  })
  #remove local files
  session$onSessionEnded(function() {
    file_d=    "/"

    files_rm=grep("VennDiagram",list.files(file_d),value=TRUE)
    sapply(1:length(files_rm),function(x){
      file.remove(paste0(file_d,"/",files_rm[x]))
    })
  })
  
  #------------------------------------------------------------------------------
  ####################### Peak distribution ###############################
  #------------------------------------------------------------------------------
  # input update
  # Update choices based on conditions
  buttonClicked <- reactiveVal(FALSE)
  
  # Update button click status when the button is pressed
  observeEvent(input$update_peak_distribution, {
    buttonClicked(TRUE)
  })
  observe({
    if(is.null(samples_bed()) || is.null(pheno_name())) {
      return()
    }
    choices <- intersect(samples_bed(), pheno_name())
    updateSelectizeInput(session, "sample_sel_peak_distr", choices = choices, selected = choices[1:min(4, length(choices))])
  })
  
  # Dynamic UI for peak distribution
  output$peak_distribution_ui <- renderUI({
    if(buttonClicked()) {
      peak_distribution_list <- lapply(1:length(input$sample_sel_peak_distr), function(x) {
        peak_distribution_UI(
          headings = input$sample_sel_peak_distr[x],
          output_info1 = paste0("peak_distribution", x),
          output_info2 = paste0("peak_distributionTSS", x),
          output_download1 = paste0("data_distribution", x),
          output_download2 = paste0("data_distributionTSS", x)
        )
      })
      tagList(peak_distribution_list)
    } else {
      tagList()  # Return an empty tag list to display nothing
    }
  })
  
  # Reactive plots and data
  peakAnnot_fig1_react <- eventReactive(input$update_peak_distribution, {
    lapply(1:length(input$sample_sel_peak_distr), function(x) {
      peakAnnot_fig1(input, output, session, file_folder = file_folder(), input_sample = input$sample_sel_peak_distr[x])
    })
  })
  peakAnnot_fig2_react <- eventReactive(input$update_peak_distribution, {
    lapply(1:length(input$sample_sel_peak_distr), function(x) {
      peakAnnot_fig2(input, output, session, file_folder = file_folder(), input_sample = input$sample_sel_peak_distr[x])
    })
  })
  peakAnnot_data1_react <- eventReactive(input$update_peak_distribution, {
    lapply(1:length(input$sample_sel_peak_distr), function(x) {
      peakAnnot_data1(input, output, session, file_folder = file_folder(), input_sample = input$sample_sel_peak_distr[x])
    })
  })
  peakAnnot_data2_react <- eventReactive(input$update_peak_distribution, {
    lapply(1:length(input$sample_sel_peak_distr), function(x) {
      peakAnnot_data2(input, output, session, file_folder = file_folder(), input_sample = input$sample_sel_peak_distr[x])
    })
  })
  
  # Render and download plots/data
  observeEvent(input$update_peak_distribution, {
    if(buttonClicked()) {  #
      lapply(1:length(input$sample_sel_peak_distr), function(x) {
        output[[paste0("peak_distribution", x)]] <- renderPlot({
          peakAnnot_fig1_react()[[x]]
        })
        output[[paste0("peak_distributionTSS", x)]] <- renderPlot({
          peakAnnot_fig2_react()[[x]]
        })
        output[[paste0("data_distribution", x)]] <- downloadHandler(
          filename = function() { paste0("PeakDistribution_", input$sample_sel_peak_distr[x], ".csv") },
          content = function(file) { write.csv(peakAnnot_data1_react()[[x]], file) }
        )
        output[[paste0("data_distributionTSS", x)]] <- downloadHandler(
          filename = function() { paste0("TF-Binding_Loci_Near_TSS_Data_", input$sample_sel_peak_distr[x], ".csv") },
          content = function(file) { write.csv(peakAnnot_data2_react()[[x]], file) }
        )
      })
    }
  })
  #-----------------------------------------------------------------------------
  ########################## Peak heatmap #####################################
  #------------------------------------------------------------------------------
  # input update
  observe({
    updateSelectizeInput(session,"sample_sel_peak_heat",choices = pheno_name(), selected = c(pheno_name()[1], pheno_name()[2],pheno_name()[3],pheno_name()[4]))
  })
  output$peak_heat = renderUI({
    peak_heat_list = lapply(1:length(input$sample_sel_peak_heat),function(x){
      peak_heat_UI(headings = paste0(input$sample_sel_peak_heat[x]),output_info = paste0("peak_heat",x) )
    })
    tagList(peak_heat_list)
  })
  peak_heat= eventReactive(input$update_peak_heat,{
    lapply(1:length(input$sample_sel_peak_heat),function(x){
      peak_heat_server(input,output,session,file_folder = file_folder(),sample_sel = input$sample_sel_peak_heat[x])
    })
  })
  observeEvent(input$update_peak_heat,{
    lapply(1:length(input$sample_sel_peak_heat),function(x){
      output[[paste0("peak_heat",x)]] = renderImage({
        peak_heat()[[x]] 
      },deleteFile = FALSE)
    })
  })
  
  #-----------------------------------------------------------------------------
  ########################## Peak meta #####################################
  #------------------------------------------------------------------------------
  # input update
  observe({
    updateSelectizeInput(session,"sample_sel_peak_meta",choices = pheno_name(), selected = c(pheno_name()[1],pheno_name()[2],pheno_name()[3],pheno_name()[4]))
  })
  output$peak_meta = renderUI({
    peak_meta_list = lapply(1:length(input$sample_sel_peak_meta),function(x){
      peak_meta_UI(headings = paste0(input$sample_sel_peak_meta[x]),output_info = paste0("peak_meta",x) )
    })
    tagList(peak_meta_list)
  })
  peak_meta= eventReactive(input$update_peak_meta,{
    lapply(1:length(input$sample_sel_peak_meta),function(x){
      peak_meta_server(input,output,session,file_folder = file_folder(),sample_sel = input$sample_sel_peak_meta[x])
    })
  })
  observeEvent(input$update_peak_meta,{
    lapply(1:length(input$sample_sel_peak_meta),function(x){
      output[[paste0("peak_meta",x)]] = renderImage({
        peak_meta()[[x]] 
      },deleteFile = FALSE)
    })
  })
  #-----------------------------------------------------------------------------
  ########################## Peak tss heatmap #####################################
  #------------------------------------------------------------------------------
  # input update
  observe({
    updateSelectizeInput(session,"sample_sel_peak_TSS_heat",choices = pheno_name(), selected = c(pheno_name()[1], pheno_name()[2],pheno_name()[3],pheno_name()[4]))
  })
  output$peak_tss_heat = renderUI({
    peak_heat_tss_list = lapply(1:length(input$sample_sel_peak_TSS_heat),function(x){
      peak_heat_tss_UI(headings = paste0("Sample:",input$sample_sel_peak_TSS_heat[x]),output_info = paste0("peak_heat_tss",x) )
    })
    tagList(peak_heat_tss_list)
  })
  peak_tss_heat= eventReactive(input$update_peak_tss_heat,{
    lapply(1:length(input$sample_sel_peak_TSS_heat),function(x){
      peak_heat_tss_server(input,output,session,file_folder = file_folder(),sample_sel = input$sample_sel_peak_TSS_heat[x])
    })
  })
  observeEvent(input$update_peak_tss_heat,{
    lapply(1:length(input$sample_sel_peak_TSS_heat),function(x){
      output[[paste0("peak_heat_tss",x)]] = renderImage({
        peak_tss_heat()[[x]] 
      },deleteFile = FALSE)
    })
  }) 
  
  #------------------------------------------------------------------------------
  ########################## peak TSS meta #####################################
  #------------------------------------------------------------------------------
  # input update
  observe({
    updateSelectizeInput(session,"sample_sel_peak_TSS_meta",choices = pheno_name(), selected = c(pheno_name()[1],pheno_name()[2],pheno_name()[3],pheno_name()[4]))
  })
  output$peak_tss_meta = renderUI({
    peak_tss_meta_list = lapply(1:length(input$sample_sel_peak_TSS_meta),function(x){
      peak_tss_meta_UI(headings = paste0("Sample:",input$sample_sel_peak_TSS_meta[x]),output_info = paste0("peak_tss_meta",x) )
    })
    tagList(peak_tss_meta_list)
  })
  peak_tss_meta= eventReactive(input$update_peak_tss_meta,{
    lapply(1:length(input$sample_sel_peak_TSS_meta),function(x){
      peak_tss_meta_server(input,output,session,file_folder = file_folder(),sample_sel = input$sample_sel_peak_meta[x])
    })
  })
  observeEvent(input$update_peak_tss_meta,{
    lapply(1:length(input$sample_sel_peak_TSS_meta),function(x){
      output[[paste0("peak_tss_meta",x)]] = renderImage({
        peak_tss_meta()[[x]] 
      },deleteFile = FALSE)
    })
  })
  

  #------------------------------------------------------------------------------------
  ######################### DEG table ##################################
  #------------------------------------------------------------------------------------
  ########################### DEG table#################################################
  # if end is smaller than start, then error
  validate_range <- reactive({
    req(input$start_deg, input$end_deg)
    
    if (input$end_deg != 0 && input$end_deg < input$start_deg) {
      return("End should be greater than start, unless it's 0.")
    }
  })
  
  observe({
    validation_result <- validate_range()
    
    if (!is.null(validation_result)) {
      showNotification(
        validation_result,
        type = "error",
        duration = 55
      )
    }
  })
  

  counts = reactive({
    dt = read_fst(paste0(file_folder(),"bw_read_counts.fst"))
    if(input$start_deg==0&input$end_deg==0&input$region_deg=="All"){
      dt = dt
    }else if(input$start_deg!=0&input$end_deg==0&input$region_deg=="All"){
      dt = subset(dt,dt$Start>=input$start_deg&dt$Chr==input$region_deg)
    } else if(input$start_deg==0&input$end_deg!=0&input$region_deg=="All"){
      dt = subset(dt,dt$End<=input$end_deg&dt$Chr==input$region_deg)
    }else{
      dt = subset(dt,dt$Start>=input$start_deg&dt$End<=input$end_deg&dt$Chr==input$region_deg)
      
    }
    dt
  })
  
  observe({
    updateSelectizeInput(session, "deg_group", choices = limma_choices(phenotype_to_use()), selected = limma_choices(phenotype_to_use())[5])
  })
  
  # test = eventReactive(input$update_DEG,{
  #   counts()
  # })
  #
  # output$deg_table = renderDT({
  #   datatable( test(),escape=F,selection = 'single')
  # })

  deg_table = eventReactive(input$update_DEG,{
    tryCatch({
      limma_table(input,output,session,counts(),pheno = phenotype_to_use(),group_info="deg_group")
    },error = function(e){
      data.frame("Need at least two genes to fit a mean-variance trend")
    })

  })
  # include geneNames and link
  deg_table_edit = eventReactive(input$update_DEG,{
   # tryCatch({
      res = deg_table()
      res$Link=paste0("<a href=","https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=", res$GeneID," target='_blank'>link</a>")
      res = data.frame(res[,c(2,1,12,3:11)])
      names(res)[which(names(res)=="log2FoldChange")] = "log2fc"
      names(res)[which(names(res)=="log.odds.of.differential.expression")] = "log.odds.of.DE"
    # },error = function(e){
    #   res = data.frame("Need at least two genes to fit a mean-variance trend")
    # })


    res
  })
  output$deg_table = renderDT({
    datatable( deg_table_edit(),escape=F,selection = 'single')

  })

  output$download_DEG =  downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("DEG_res", ".csv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write.csv(deg_table(), file)
    }
  )
  
  #------------------------------------------------------------------------------------
  ###################### region plot (exclude from the app)############################
  #------------------------------------------------------------------------------------
  observe({
    trigger = input$mainnav_region
    updateSelectizeInput(session,"gene_region",choices = c("",counts()$GeneName),selected="",server=T)
  })
  # grayout part of the option
  observeEvent(input$gene_region, {
    if(input$gene_region != "") {
      disable('region')
      disable("start")
      disable("end")
    } else {
      enable('region')
      enable("start")
      enable("end")
    }
  }, ignoreNULL = T)
  #------------------------------------------------------------------------------------
  ######################### DEG figure ##################################
  #------------------------------------------------------------------------------------
  geneID = reactive({
    res = deg_table()
    paste0(res$GeneName,"_",res$GeneID)
  })
  observe({
    updateSelectizeInput(session,"deg_fill",choices  = geneID(),server=T)
    
  })
  deg_table_fig = reactive({
    res = deg_table()
#    res$gene_name =paste0(res$GeneName,"_",res$Chr,"_",res$Start,"_",res$End,"_",res$GeneID)
    res$gene_name =paste0(res$GeneName,"_",res$GeneID)
    
    rownames(res) = res$gene_name
    res = data.frame(res[,6:(ncol(res)-1)])
    for(i in 1:6){
      res[,i] = as.numeric(as.character(res[,i]))
    }
    res
  })
  
  deg_fig = eventReactive(input$update_MA,{
    deg_figServer(input,output,session,diff_result=deg_table_fig(),folder_direct=folder_direct(),fill_ma="deg_fill",title = input[["deg_group"]])
  })
  output$deg_fig = renderPlotly({
    deg_fig()
  })
  
  #---------------------------------------------------------------------------------
  ##################### DEG region plot ##################################
  #---------------------------------------------------------------------------------
  observe({
   # updateSelectizeInput(session,"sel_sample_regionDeg",choices = samples(),selected = c(samples()[3],samples()[2]),server =T)
    group_info1 = str_split(input[["deg_group"]],":")[[1]][1]
    comp_group = str_split(input[["deg_group"]],":")[[1]][2]
    comp_group_sub = str_split(comp_group," vs ")[[1]]
    sample1 = phenotype_to_use()$X[phenotype_to_use()[[group_info1]]==comp_group_sub[1]]
    sample2 = phenotype_to_use()$X[phenotype_to_use()[[group_info1]]==comp_group_sub[2]]
    updateSelectizeInput(session,"sel_sample_regionDeg",choices = c(sample1,sample2),selected = c(sample1[1],sample2[1]),server =T)
  })
  
  observe({
    trigger = input$mainnav_deg
    counts = counts()
    updateSelectizeInput(session,"gene_region_deg",choices = c("",counts$GeneName[!is.na(counts$GeneName)]),selected="",server=T)
  })
  
  observeEvent(input$gene_region_deg, {
    if(input$gene_region_deg != "") {
      disable('chr_region_deg')
      disable("start_deg_region")
      disable("end_deg_region")
    } else {
      enable('chr_region_deg')
      enable("start_deg_region")
      enable("end_deg_region")
    }
  }, ignoreNULL = T)
  
  pheno = reactive({
    dt = fread(paste0(file_folder(),"phenotype.csv"))
    dt = data.frame(dt)
    dt
  })
  
  sample1 = reactive({
    group_info1 = str_split(input[["deg_group"]],":")[[1]][1]
    comp_group = str_split(input[["deg_group"]],":")[[1]][2]
    comp_group_sub = str_split(comp_group," vs ")[[1]]
    
    sample1 = pheno()$X[pheno()[[group_info1]]==comp_group_sub[1]]
    sample1
  })
  sample2 = reactive({
    group_info1 = str_split(input[["deg_group"]],":")[[1]][1]
    comp_group = str_split(input[["deg_group"]],":")[[1]][2]
    comp_group_sub = str_split(comp_group," vs ")[[1]]
    
    sample2 = pheno()$X[pheno()[[group_info1]]==comp_group_sub[2]]
    sample2
  })
  
  
  deg_coverage = eventReactive(input$update_region_deg,{
    # chr = input$chr_region_deg;start = input$start_deg_region;end=input$end_deg_region
    # s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
    
    
    if(input$gene_region_deg==""){
      chr = input$chr_region_deg;start = input$start_deg_region;end=input$end_deg_region
      s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
    }else{
      counts = counts()
      chr = counts$Chr[which(counts$GeneName==input$gene_region_deg)]
      start = counts$Start[which(counts$GeneName==input$gene_region_deg)]
      end = counts$End[which(counts$GeneName==input$gene_region_deg)]
      s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
      
    }
    
    setwd(paste0(file_folder(),"bigwig/"))
    # Load the ChIP-seq data for the samples
    normal_samples = paste0(sample1(),".bw")
    tumor_samples = paste0(sample2(),".bw")
    
    
    # Create a list to store coverage data for each sample
    coverage_list <- list()
    
    # Load and preprocess data for normal samples
    for (bw_file in normal_samples) {
      bw <- GRanges(import.bw(bw_file, selection = BigWigSelection(s)))
      bw <- coverage(bw, weight = 'score')
      y <- as.vector(unlist(bw[[seqnames(s)]][start(s):end(s)]))
      y <- runmean(y, 30, endrule = 'keep')
      coverage_list[[bw_file]] <- y
    }
    
    # Load and preprocess data for tumor samples
    for (bw_file in tumor_samples) {
      bw <- GRanges(import.bw(bw_file, selection = BigWigSelection(s)))
      bw <- coverage(bw, weight = 'score')
      y <- as.vector(unlist(bw[[seqnames(s)]][start(s):end(s)]))
      y <- runmean(y, 20, endrule = 'keep')
      coverage_list[[bw_file]] <- y
    }
    
    # Create colors for normal and tumor samples
    colors <- c(rep(brewer.pal(3, 'Set1')[1], length(normal_samples)),
                rep(brewer.pal(3, 'Set1')[2], length(tumor_samples)))
    
    # Combine coverage data for all samples
    combined_data <- do.call(cbind, coverage_list)
    
    combined_df <- as.data.frame(combined_data)
    combined_df$Position <- seq(start, end, length.out = nrow(combined_df))
    
    melted_df <- melt(combined_df, id.vars = "Position", variable.name = "Sample", value.name = "Coverage")
    # Filter out zero values
    melted_df <- melted_df[melted_df$Coverage != 0, ]
    # Add phenotype group information
    melted_df$Group <- ifelse(melted_df$Sample %in% normal_samples, 'group1', 'group2')
    
    # Plot using ggplot2
    p = ggplot(melted_df, aes(x = Position, y = Coverage, group = Sample, color = Sample, linetype = Group)) +
      geom_line(linewidth=1) +
      theme_minimal() +
      labs(x = seqnames(s), y = "ChIP-seq coverage") +
      scale_color_brewer(palette = "Set1") +
      theme(legend.title = element_blank(), legend.position = "bottom") +
      guides(col = guide_legend("Phenotype Group",nrow=4,byrow=TRUE), linetype = guide_legend("Sample"))
    
    p
  })
  
  
  
  output$region_plot_deg = renderPlot({
    deg_coverage()
  })
  
  #-------------------------------------------------------------------------------
  #################### DEG region for each sample #####################
  #------------------------------------------------------------------------------
  deg_coverage_each_fig = eventReactive(input$update_region_deg,{
    
    group_info1 = str_split(input[["deg_group"]],":")[[1]][1]
    comp_group = str_split(input[["deg_group"]],":")[[1]][2]
    comp_group_sub = str_split(comp_group," vs ")[[1]]
    
    if(input$gene_region_deg==""){
      chr = input$chr_region_deg;start = input$start_deg_region;end=input$end_deg_region
      s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
    }else{
      counts = counts()
      chr = counts$Chr[which(counts$GeneName==input$gene_region_deg)]
      start = counts$Start[which(counts$GeneName==input$gene_region_deg)]
      end = counts$End[which(counts$GeneName==input$gene_region_deg)]
      s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
      
    }
    
    setwd(paste0(file_folder(),"bigwig/"))
    # Load the ChIP-seq data for the samples
    normal_samples = paste0(sample1(),".bw")
    tumor_samples = paste0(sample2(),".bw")
    
    
    
    
    
    # Create a list to store coverage data for each sample
    coverage_list <- list()
    
    # Load and preprocess data for normal samples
    for (bw_file in normal_samples) {
      bw <- GRanges(import.bw(bw_file, selection = BigWigSelection(s)))
      bw <- coverage(bw, weight = 'score')
      y <- as.vector(unlist(bw[[seqnames(s)]][start(s):end(s)]))
      y <- runmean(y, 30, endrule = 'keep')
      coverage_list[[bw_file]] <- y
    }
    
    # Load and preprocess data for tumor samples
    for (bw_file in tumor_samples) {
      bw <- GRanges(import.bw(bw_file, selection = BigWigSelection(s)))
      bw <- coverage(bw, weight = 'score')
      y <- as.vector(unlist(bw[[seqnames(s)]][start(s):end(s)]))
      y <- runmean(y, 20, endrule = 'keep')
      coverage_list[[bw_file]] <- y
    }
    
    
    
    colors <- c(rep(brewer.pal(3, 'Set1')[1], length(normal_samples)),
                rep(brewer.pal(3, 'Set1')[2], length(tumor_samples)))
    
    names(coverage_list) <- sub(".bw","",names(coverage_list))
    
    normal_samples =  sub(".bw","",normal_samples)
    tumor_samples =  sub(".bw","",tumor_samples)
    
    sample_groups <- ifelse(names(coverage_list) %in% normal_samples, comp_group_sub[[1]], comp_group_sub[[2]])
    names(sample_groups) = names(coverage_list)
    
    # Calculate global y-limits
    global_min <- min(sapply(coverage_list, min))
    global_max <- max(sapply(coverage_list, max))
    global_ylims <- c(global_min, global_max)
    
    # Define the layout for plots and legend
    num_samples <- length(input$sel_sample_regionDeg)
    par(mfrow = c(num_samples + 1, 1), mar = c(2, 2, 2, 2))
    
    # Loop through each sample and create a plot
    for (i in input$sel_sample_regionDeg) {
      sample_data <- coverage_list[[i]]
      # Generate a sequence of genomic positions
      genomic_positions <- seq(from = start, to = end, length.out = length(sample_data))
      
      # Set the margin for each plot
      par(mar = c(2, 2, 2, 2))
      plot(genomic_positions, sample_data, type = 'l',
           col = ifelse(sample_groups[i] == comp_group_sub[[1]], 'blue', 'red'),
           main = names(coverage_list)[i], xlab = "Genomic Position", ylab = "Coverage",
           ylim = global_ylims, xlim = c(start, end), xaxt = "n")
      
      
      # Add custom x-axis ticks
      if (i == tail(input$sel_sample_regionDeg, n = 1)) {
        axis(1, at = seq(from = start, to = end, length.out = 5),
             labels = formatC(seq(from = start, to = end, length.out = 5), format = "d", digits = 0))
      }
      
    }
    
    # Add the legend in its designated area
    plot.new()
    legend("center", legend = c(comp_group_sub[[1]], comp_group_sub[[2]]), col = c("blue", "red"), lty = 1, cex = 1, horiz = TRUE)
    
    
  })
  
  output$region_plot_DEG_each = renderPlot({
    deg_coverage_each_fig()
  })
  
  #---------------------------------------------------------------------------------
  ############# overlaped samples in region plot ############
  #----------------------------------------------------------------------------------
  deg_coverage_overlap = eventReactive(input$update_region_deg,{
    # chr = input$chr_region_deg;start = input$start_deg_region;end=input$end_deg_region
    # s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
    
    
    if(input$gene_region_deg==""){
      chr = input$chr_region_deg;start = input$start_deg_region;end=input$end_deg_region
      s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
    }else{
      counts = counts()
      chr = counts$Chr[which(counts$GeneName==input$gene_region_deg)]
      start = counts$Start[which(counts$GeneName==input$gene_region_deg)]
      end = counts$End[which(counts$GeneName==input$gene_region_deg)]
      s = GRanges(chr,IRanges(as.numeric(start),as.numeric(end)))
      
    }
    
    setwd(paste0(file_folder(),"bigwig/"))
    # Load the ChIP-seq data for the samples
    normal_samples = paste0(sample1(),".bw")
    tumor_samples = paste0(sample2(),".bw")
    
    
    # Create a list to store coverage data for each sample
    coverage_list <- list()
    
    # Load and preprocess data for normal samples
    for (bw_file in normal_samples) {
      bw <- GRanges(import.bw(bw_file, selection = BigWigSelection(s)))
      bw <- coverage(bw, weight = 'score')
      y <- as.vector(unlist(bw[[seqnames(s)]][start(s):end(s)]))
      y <- runmean(y, 30, endrule = 'keep')
      coverage_list[[bw_file]] <- y
    }
    
    # Load and preprocess data for tumor samples
    for (bw_file in tumor_samples) {
      bw <- GRanges(import.bw(bw_file, selection = BigWigSelection(s)))
      bw <- coverage(bw, weight = 'score')
      y <- as.vector(unlist(bw[[seqnames(s)]][start(s):end(s)]))
      y <- runmean(y, 20, endrule = 'keep')
      coverage_list[[bw_file]] <- y
    }
    
    colors <- c(rep(brewer.pal(3, 'Set1')[1], 1),
                rep(brewer.pal(3, 'Set1')[2], 1))
    # Combine coverage data for all samples
    combined_data <- do.call(cbind, coverage_list[paste0(input$sel_sample_regionDeg,".bw")])
    
    # Create the plot with different colors for each sample
    matplot(start(s):end(s), combined_data, type = 'l', ylim = c(0, max(combined_data) * 2),
            col = c(colors), lwd = 2, xlab = seqnames(s), ylab = 'ChIP-seq coverage')
    
    # create region plot for each sample indicidually
    
    # Add a legend to differentiate normal and tumor samples
    legend('topright', legend = c(input$sel_sample_regionDeg), col = unique(colors), lwd = 2, bty = 'n')
    
  })
  
  output$region_plot_DEG_overlap = renderPlot({
    deg_coverage_overlap()
  })  
}



# Run the application 
shinyApp(ui = ui, server = server)