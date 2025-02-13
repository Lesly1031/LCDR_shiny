#===============================================================================
# Create ShinyApp with enable uploading data
#===============================================================================
packages = c("shiny","ggplot2","ggpubr","tidyverse","DT","data.table","plotly","shinyjs","shinyBS","pheatmap","DESeq2","edgeR",
             "msigdbr","shinycustomloader","jose","EBSeq","webshot","shinyWidgets","org.Mm.eg.db","grid","xCell","immunedeconv",
             "org.Hs.eg.db","ggrepel","fst","data.table","DEGreport","nipals","enrichR","stringi","rhandsontable")


#install.packages("rhandsontable")

#install.packages("rhandsontable")

#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    # some packages must be installed in seperately
    if (x %in% c("ComplexHeatmap","DESeq2","edgeR","clusterProfiler","enrichplot","pathview","EBSeq")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(x)
      
    } 
    #else if(x=="shinyTree"){
    # devtools::install_github("shinyTree/shinyTree")
    #}
    else {
      install.packages(x, dependencies = TRUE,repos = "http://cran.us.r-project.org")
    }
  }
  library(x, character.only = TRUE)
})

source("/Shiny_uploadData/sample_table.R")
source("/pca_plot.R")
source("/DEG_table.R")
source("/DEG_fig.R")
source("/heat_deg.R")
source("/heat_kegg.R")
source("/box_plot.R")
source("/scatter_plot.R")



# Define UI for application that draws a histogram
ui <- fluidPage(
  useShinyjs(),
  tags$script("
    Shiny.addCustomMessageHandler('resetValue', function(variableName) {
      Shiny.onInputChange(variableName, null);
    });
  "),
  theme = "style/style.css",
  tagList(
    div(id="header",
        
        h4("  LCDR - Visualization"),
        
    ),
    tags$img(align="right",src="images/logo_Moffitt.jpg",height="45px"),
    tabsetPanel(type="tabs",id="mainnav", 
                tabPanel(
                  "Upload Data",
                  fluidRow(
                    column(2,
                           selectizeInput("study_type","Select the study type",c("RNAseq","Metabolomics","Proteomics","lipidomics"),"Metabolomics"),
                             div(fileInput('file_omics',
                                           span(
                                             list(HTML("<p><abbr title='Please upload a CSV file where the first column contains gene names labeled as 'geneName,' and the remaining columns represent samples'>Upload your phenotype file (csv file)....</abbr></p>"))
                                           ),
                                           accept = c(
                                             '.csv'
                                           )),
                                 style="font-size:100%;"
                             ),
                             div(fileInput('file_pheno',
                                           span(
                                             list(HTML("<p><abbr title='Please upload your phenotype file in CSV format. The first column should be labeled 'SampleID' and contain sample identifiers, with the subsequent columns containing phenotype information'>Upload your phenotype file (csv file)....</abbr></p>"))
                                           ),
                                           accept = c(
                                             '.csv'
                                           )),
                                 style="font-size:100%;"
                             ),
                             actionButton("load_omics", "View omics Example"),
                             bsModal("omics_example_window", "", "load_omics", size = "large",
                                     withLoader(
                                       
                                       DTOutput("omics_example"),type="html",loader="loader1")),
                             
                             div(style = "padding-top: 20px;"),
                             actionButton("load_pheno", "View phenotype Example"),
                             bsModal("pheno_example_window", "", "load_pheno", size = "large",
                                     withLoader(DTOutput("pheno_example"),type="html",loader="loader1"))
                           
                           ),
                    column(8,
                           div(id="box1",
                               panel(
                                 status="info",heading = "View the uploaded omics data",
                                 div(style='max-width: 100%; height: 500px; width: auto;overflow-y: scroll;',
                                     withLoader(DTOutput("upload_omics_view"),type="html",loader="loader1")
                                 )
                               )  
                           )
                           )
                  )
                ),
                tabPanel("Phenotype",
                         fluidRow(
                           column(2,
                                  actionButton("instruction_pheno","Tab instructions")
                           ),
                           column(10,
                                  uiOutput("all_sample"),
                                  DTOutput('test_sample')
                           )
                         )
                         
                ),
                tabPanel("PCA",
                         fluidRow(
                           column(2,
                                  actionButton("update_pca","Run")
                                  
                           ),
                           
                           column(10,
                                  uiOutput("all_pca_results")
                           )
                           
                         )
                ),
                
                tabPanel("DEG",
                         tabsetPanel(type="tabs",id="mainnav_sub",
                                     tabPanel("DEG test result",
                                              fluidRow(
                                                column(2,
                                                       br(),
                                                       p("Click the button below after selecting the comparison group from the drop down gear button"),
                                                       
                                                       actionButton("update_DEG","Run")
                                                      # br(),
                                                      # p('If you select an RNAseq study, you will have the option to download the normalized data following voom transformation. Simply click on the gear icon in the top left corner of the study panel and select "Download voom data." Additionally, by clicking on "Download DEG table," you can download the table of differentially expressed genes (DEG).')
                                                ),
                                                column(10,
                                                       uiOutput("all_deg_results")
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
                                   
                                     tabPanel("Heatmap for DEG",
                                              fluidRow(
                                                column(2,
                                                       numericInput("top_n_heat_deg","Show top n features",value=25),
                                                       checkboxInput("show_row_names_deg",
                                                                     "Show omics and sample names", value = FALSE),
                                                       # checkboxInput("show_column_names_deg",
                                                       #               "Show sample names", value = FALSE),
                                                       actionButton("update_heat_deg","Run")
                                                ),
                                                column(10,uiOutput("all_heat_deg"),
                                                       DTOutput("test_degheat")
                                                       
                                                       )
                                              )
                                     ),
                                     tabPanel("KEGG for DEG",
                                              fluidRow(
                                                column(2,
                                                       actionButton("instruction_kegg","Tab instructions"),
                                                       # selectInput("which_p_kegg","Choose p value or adjusted p",choices = c("pval","adjust.p"),selected="pval"),
                                                       # numericInput("p_threshold_kegg", "pval/adj.pval threshold", 0.05,min=0.0001,max=0.1,step=0.05),
                                                       # numericInput("fc_threshold_kegg", "Select fold change threshold", value=1.5),
                                                       checkboxInput("show_row_names_kegg",
                                                                     "Show gene names", value = FALSE),
                                                       #selectizeInput("omics_filter_kegg","Filter number of genes/metanbolites/proteomics",c("Yes","No"),"No"),
                                                       # conditionalPanel(
                                                       #   condition="input.omics_filter_kegg=='Yes'",
                                                       #   numericInput("number_omics_kegg","Only showing pathways with at least the selected number of genes/metabolites/proteomics",value=20,min=3,step=1)
                                                       # ),
                                                       p("Click the button below after selecting the pathway from the drop down gear button"),
                                                       actionButton("update_heat_kegg","Run")
                                                ),
                                                column(10,uiOutput("all_keggheat_results"))
                                               # DTOutput("test")
                                              )
                                     )
    
                         )
                         
                ),
                tabPanel("Box plot",
                         fluidRow(
                           column(2,
                                  #radioButtons("box_meta","Only include identified metabolites? (only for metabolomics studies)",c("Yes","No"),"No"),
                                  actionButton("update_box","Run")
                                  
                           ),
                           column(10,
                                  uiOutput("all_boxplot_results")
                           )
                         )
                ),
                
                # tabPanel("Heatmap",
                #          tabsetPanel(type="tabs",id="mainnav_kegg_any",
                #            tabPanel("Heatmap",
                #                     fluidRow(
                #                       column(2,
                #                              
                #                              # checkboxInput("cluster_rows",
                #                              #               "Cluster rows", value = TRUE),
                #                              # checkboxInput("cluster_columns",
                #                              #               "Cluster columns", value = TRUE),
                #                              # radioButtons("heat_meta","Only include identified metabolites? (only for metabolomics studies)",c("Yes","No"),"No"),
                #                              
                #                              checkboxInput("show_row_names",
                #                                            "Show omics names", value = FALSE),
                #                              # checkboxInput("show_column_names",
                #                              #               "Show sample names", value = FALSE),
                #                              
                #                              p("Shows the top 25 markers from DEG results by default"),
                #                              actionButton("update_heat","Run")
                #                              
                #                              # checkboxInput("z_score",
                #                              #               "Z-score transformation", value = TRUE)
                #                       ),
                #                       column(10,uiOutput("all_heat_results"))
                #                     )
                #            ),
                #            tabPanel("KEGG for selected genes",
                #                     fluidRow(
                #                       column(2,
                # 
                #                              checkboxInput("show_row_names_kegg_any",
                #                                            "Show gene names", value = FALSE),
                #                            
                #                              p("Click the button below after selecting the pathway from the drop down gear button"),
                #                              actionButton("update_heat_kegg_any","Run")
                #                       ),
                #                       column(10,uiOutput("all_keggheat_results_any"))
                #                     )
                #            )
                # 
                #            
                #            
                #          )
                #          
                #          
                #          
                #          
                # ),
                tabPanel("Scatterplot",
                         tabsetPanel(type="tabs",
                                     tabPanel("For each study",
                                              fluidRow(
                                                column(2,
                                                       
                                                       #  radioButtons("scatter_meta","Only include identified metabolites? (only for metabolomics studies)",c("Yes","No"),"No")
                                                       actionButton("update_scatter","Run")
                                                       
                                                ),
                                                column(10,
                                                       uiOutput("all_scatter_results"),
                                                       textOutput("text")
                                                       #dataTableOutput("test")
                                                )
                                                
                                              )
                                     )
                         
                         )
                         
                ),
                
                # tabPanel("Immune decovolution",
                #          tabsetPanel(type="tabs",id="mainnav_immunedeconv", 
                #            tabPanel("Immune deconv table",
                #                     fluidRow(
                #                       column(2,
                #                              p("This tab displays heatmaps for immune deconvolution in RNA-seq studies. The data has been transformed using z-scores"),
                #                              br(),
                #                              p("Click 'Run' button to get the table after selection"),
                #                              actionButton("update_table_immuneDeconv","Run"),
                #                              
                #                       ),
                #                       column(10,
                #                              uiOutput("all_immuneDeconv_table"),
                #                              DTOutput("test_deconv")
                #                       )
                #                     )
                #                     ),
                #            tabPanel("Immune deconv heatmap",
                #                     fluidRow(
                #                       column(2,
                #                              br(),
                #                              p("Click 'Run' button to get the heatmap after selection"),
                #                              actionButton("update_heat_immuneDeconv","Run"),
                #                              
                #                       ),
                #                       column(10,
                #                              uiOutput("all_immuneDeconv_heat")
                #                       )
                #                     )
                #                     ),
                #            tabPanel("immune deconv boxplot",
                #                     fluidRow(
                #                       column(2,
                #                              br(),
                #                              p("Click 'Run' button to get the boxplot after selection"),
                #                              actionButton("update_box_immuneDeconv","Run"),
                #                              
                #                       ),
                #                       column(10,
                #                              uiOutput("all_immuneDeconv_box")
                #                       )
                #                     )
                #                     )
                #          )
                # 
                #          )

    )
  )
  
  
  
)

options(shiny.maxRequestSize=30*1024^2) 


server <- function(input, output,session) { 
  
  folder_direct <- reactive({
    file = "/Shiny_uploadData/"
  })
  
#------------------------------------------------------------------------------
# upload dataset #
#------------------------------------------------------------------------------
  omics_example = reactive({
    fread(paste0(folder_direct(),"www/example_data/omics_example.csv"))
    
  })
  pheno_example = reactive({
    fread(paste0(folder_direct(),"www/example_data/phenotype_example.csv"))
  })
  
  omics <- reactive({
    infile <- input$file_omics
    if (is.null(infile)) {
      omics = fread(paste0(folder_direct(),"www/example_data/omics_example.csv"))
    }else{
      omics= fread(infile$datapath)
    }
    omics = data.frame(omics,check.names=F)
    colnames(omics)[1] = "geneName"
    return(omics)
  })
  
  pheno = reactive({
    infile <- input$file_pheno
    if (is.null(infile)) {
      pheno = fread(paste0(folder_direct(),"www/example_data/phenotype_example.csv"))
    } else{
      pheno = fread(infile$datapath)
    }
    pheno = data.frame(pheno,check.names=F)
    colnames(pheno)[1] = "SampleID"
    return(pheno)
  })

  output$upload_omics_view = renderDT({
    data.table(omics(),options = list(scrollX = TRUE
    ))
    
    
  })
 
  ##################### click the button to show the example omics (two small dataset)#####################
  output$omics_example = renderDT({
    datatable(omics_example(),options = list(scrollX = TRUE))
  })
  
  output$pheno_example = renderDT({
    pheno_example()
  })
  
  #===============================================================================
  ########################## phenotype data file ##############################
  #===============================================================================
  output$all_sample<-renderUI({
    sample_UI(output_info = paste0("sample","1"))
  })
  
  
  observe({
      output[[paste0("sample","1")]]<-renderRHandsontable({
        sample_server(input,output,session,pheno_data = pheno())
      })
    
  })
  
  phenotype_to_use <- reactive({
      dt_class=hot_to_r(input[[paste0("sample","1")]])
      dt_class=data.frame(dt_class,check.names = F)
      dt_class = dt_class[complete.cases(dt_class$SampleID),]
      names(dt_class)[1] = "X"
      dt_class
    
  })
  
  genotype <- reactive({
      dt = omics()
      names(dt)[1] = "ID"
      dt
  })

  geneID_ID = reactive({
    genotype()$ID
  })
  
  study_model = reactive({
    lapply(1:length(folder_direct()),function(x){
      id = genotype()[[x]][["ID"]]
      id = sapply( strsplit(id, "_"), "[", 1)
      uppercase_id <- all(grepl("^[[:upper:]]+$", id[2]))
      if(isTRUE("all_uppercase")){
        study_model="human"
      }else{
        study_model=  "mouse"
      }
      study_model
    })
  })
  
  
  # 
  
  onclick('instruction_pheno', showModal(modalDialog(size="l",fade=TRUE,
                                                     title = "",
                                                     renderUI(
                                                       tags$img(src="images/pheno_info.jpg",width="800px")
                                                     ),easyClose = TRUE
  ))) 
  
  #==============================================================================
  ############################# PCA ############################################
  #==============================================================================
  
  ############################PCA#############################
  output$all_pca_results<-renderUI({
      # genotype = values[[paste0("genotype",x)]]
      pca_UI(headings = paste0(input$study_type," - PCA"),output_info = paste0("pca",1),group_info1 = paste0("group_pca1",1),group_info2 = paste0("group_pca2",1),download_info=paste0("downloadpca",1),
             phenotype = phenotype_to_use(),study_type=input$study_type)
    
  })
  
  # ###########pca server###########
  pca_res = eventReactive(input$update_pca,{

      pca_server(input,output,session,data=omics(),data_class=phenotype_to_use(),
                 group_info1=paste0("group_pca1","1"),group_info2=paste0("group_pca2","1"),study_type=input$study_type)

  })

  observe({
      output[[paste0("pca","1")]]<-renderPlotly({
        p =  pca_res()
        ggplotly(p, tooltip="text") %>%
          layout(legend = list(orientation = 'h', x = 1.02, y = 1))
      })

  })
  # 
  # Download figure
  observe({
      output[[paste0("downloadpca",1)]]<-downloadHandler(
        filename = function() {
          paste0("PCA_",input$study_type,'.png', sep = '')

        },
        content = function(file){
          ggsave(file, pca_res(),type="cairo-png",width = 15, height = 15,units="cm")
        }
      )
  })
      ########################### DEG table#################################################
      ########################### DEG table#################################################
      output$all_deg_results<-renderUI({
          # genotype = values[[paste0("genotype",x)]]
          limma_UI(headings=input$study_type,output_info = paste0("limma",1),
                   group_info = paste0("group_limma",1),
                   download_info=paste0("downloadlimma",1),
                   phenotype = phenotype_to_use(), 
                   study_type = input$study_type)
        
      })
  
  ######## DEG table server ##########
  deg_res <- eventReactive(input$update_DEG,{
      res = limma_server(input,output,session,data=genotype(),data_class = phenotype_to_use(),
                   group_info = paste0("group_limma",1),study_type =input$study_type)
      deg_res = res
      deg_res$rowID = rownames(deg_res)
      rownames(deg_res) = NULL
      deg_res$gene_name = deg_res$rowID
      
      deg_res = data.frame(deg_res[,c("gene_name","log2FoldChange","baseMean","pvalue","padj")],check.names=F)
      
      deg_res
  })
  observe({
      output[[paste0("limma",1)]]<-renderDT({
        deg_res = deg_res()
        # deg_res$rowID =  sapply( str_split(rownames(deg_res), "_"), "[", 1)
 
        
        # #deg_res = data.frame(deg_res[,c("rowID",names(deg_res)[-ncol(deg_res)]  )])
        # 
        # #        if(grepl("metabolomics",folder_direct()[[x]])|grepl("lipidomics",folder_direct()[[x]])){
        # if(grepl("metabolomics",folder_direct()[[x]])){
        #   
        #   # names = sapply(strsplit(deg_res$rowID,"_"),"[",1) 
        #   # IDs = unlist(lapply(1:length(deg_res$rowID),function(x){
        #   #   paste(strsplit(deg_res$rowID[x],"_")[[1]][2:3],collapse = "_")
        #   # }))
        #   # deg_res$Names= ifelse(names=="X",deg_res$rowID,names)
        #   # deg_res$IDs= ifelse(names=="X",deg_res$rowID,IDs)
        #   split_sample_names <- function(names) {
        #     # The first part up to the polarity (pos or neg)
        #     part1 <- str_extract(names, ".*(?=_pos|_neg)")
        #     # The second part starting from the polarity (pos or neg)
        #     part2 <- str_extract(names, "(pos|neg)_.+")
        #     # Combine the parts into a data frame
        #     data.frame(part1, part2)
        #   }
        #   
        #   # Apply the function to split the names
        #   split_names_df <- split_sample_names(deg_res$rowID)
        #   names = split_names_df$part1
        #   IDs = split_names_df$part2
        #   deg_res$Names= ifelse(names=="X",deg_res$rowID,names)
        #   deg_res$IDs= ifelse(names=="X",deg_res$rowID,IDs)
        #   
        # }else if(grepl("lipidomics",folder_direct()[[x]])){
        #   deg_res$Names = deg_res$rowID
        #   deg_res$IDs = ""
        # }else{
        #   deg_res$Names = sapply(strsplit(deg_res$rowID,"_"),"[",1)
        #   deg_res$IDs = sapply(strsplit(deg_res$rowID,"_"),"[",-1)
        #   
        # }
        #deg_res$Names = deg_res$rowID
        
        #deg_res = data.frame(deg_res[,c("Names","log2FoldChange","baseMean","pvalue","padj")],check.names=F)
        datatable(
          deg_res ,
          rownames = FALSE,plugins = "ellipsis",
          extensions = 'Buttons',
          options = list(
            autoWidth = FALSE, scrollX = TRUE,
            columnDefs = list(list(
              width = "10px", targets = "_all",
              render = JS("$.fn.dataTable.render.ellipsis( 30, false )")
            )))
        )%>%
          formatSignif(columns = c('pvalue', 'padj'), digits = 3)
        
      })
    })
  
  ########################### download DEG tables ################################
  # Download TABLE 
  observe({
      output[[paste0("downloadlimma",1)]]<-downloadHandler(
        filename = function() {
          paste0("DEG_",input$study,'.csv', sep = '')
          
        },
        content = function(file){
          fwrite(deg_res(),file,row.names = T)
        }
      )
    
  })
 
######################## DEG figure #############################
  geneID = reactive({
    res = deg_res()
    paste0(res$gene_name)
  })
  
  observe({
    updateSelectizeInput(session,"deg_fill",choices  = geneID(),server=T)
    
  })
  deg_table_fig = reactive({
    res = deg_res()
    #    res$gene_name =paste0(res$GeneName,"_",res$Chr,"_",res$Start,"_",res$End,"_",res$GeneID)

    for(i in 2:ncol(res)){
      res[,i] = as.numeric(as.character(res[,i]))
    }
    rownames(res) = res$Names
    res
  })
  
  deg_fig = eventReactive(input$update_MA,{
    deg_figServer(input,output,session,diff_result=deg_table_fig(),fill_ma="deg_fill",title = input[["deg_group"]])
  })
  output$deg_fig = renderPlotly({
    ggplotly(deg_fig(),tooltip="text")
  })
  
  ####################### Heatmap DEG ######################################
  
  output$all_heat_deg<-renderUI({
      heatmapDegUI(headings = input$study_type,output_info = paste0("heat_deg",1),color_annot=paste0("color_annot_heat_deg",1))

  })
  
  
  ###########heatmap server###########
  # heat_res = eventReactive(input$update_heat,{
  observeEvent(input$update_heat_deg,{

      output[[paste0("heat_deg",1)]]<-renderPlot({
        heatmapDeg_server(input,output,session,data=genotype(),data_class=phenotype_to_use(),group_info=paste0("group_limma",1),deg_res = deg_res(),group_info_deg=paste0("group_limma",1))
  
    })
  })
  #=========================================================================================
  ####################### Heatmap KEGG #################################################
  # Have to Work on KEGG for metaboliomics 
  #=========================================================================================
  # Make the tab not visible for lipdomics
  # Make the tab not visible for lipdomics
  tab_kegg_visible = reactiveVal(TRUE)
  
  observe({
      if(grepl("Lipidomics",input$study_type)|grepl("Proteomics",input$study_type)){ ################# Need to change on server
        tab_kegg_visible(FALSE)
      }else{
        tab_kegg_visible(TRUE)
      }
  })
  
  observe({
    if (tab_kegg_visible()) {
      showTab(inputId = "mainnav_sub", target = "KEGG for DEG")
    } else {
      hideTab(inputId = "mainnav_sub", target = "KEGG for DEG")
    }
  })
  
  ######################
  onclick('instruction_kegg', showModal(modalDialog(size="l",fade=TRUE,
                                                    title = "",
                                                    renderUI(
                                                      #tags$img(src="images/kegg_info.jpg",width="700px")
                                                      p("This tab provides a graphical representation of the relative abundance or expression levels of biological molecules such as gene, proteins, or metabolites across different samples or conditions. The KEGG (Kyoto Encyclopedia of Genes and Genomes) database contains information on biological pathways and molecular interactions that can aid in comprehending the functional and regulatory mechanisms of cellular processes.  KEGG heatmaps are generated by first mapping the RNAseq, proteomics or metabolomics data to KEGG pathways, and then visualizing the expression levels of the mapped genes, proteins or metabolites as a heatmap.")
                                                    ),easyClose = TRUE
  ))) 
  
  output$all_keggheat_results<-renderUI({
      kegg_heatUI(headings = input$study_type,output_info = paste0("heat_kegg",1),output_info_tb = paste0("heat_kegg_table",1),kegg_path = paste0("kegg_path",1) ,group_info = paste0("group_heat_kegg",1),
                  phenotype = phenotype_to_use(),download_info=paste0("downloadheatkegg",1),gene_pathway = paste0("gene_pathway",1),color_annot=paste0("color_annot",1)
      )
    
  })
  enrich_score<-reactive({
      enrich_res_tb(input,output,session,deg_data = data.frame(deg_res()),study_type = input$study_type)

  })
  
  
  
  observe({
      trigger = input$mainnav_sub
      enrich_res = enrich_score()
      #    if(input$omics_filter_kegg=="No"){
      enrich_res = data.frame(enrich_res)
      # }else{
      #   enrich_res = subset(enrich_res,enrich_res$length_genes<=input$number_omics_kegg&enrich_res$length_genes>=3)
      # }
      
      kegg_sel_observe(session,input_info =paste0("kegg_path",1),enrich_res = enrich_res)
    
  })
  
  
  #kegg_heat = eventReactive(input$update_heat_kegg,{
  observeEvent(input$update_heat_kegg,{
    
      enrich_data = enrich_score()
      phenotype_data = phenotype_to_use()
      omics_data = genotype()
      
      output[[paste0("heat_kegg",1)]] = renderPlot({
        kegg_heatServer(input,output,session,enrich_data =enrich_data,omics_data =omics_data ,kegg_path= paste0("kegg_path",1),group_info=paste0("group_heat_kegg",1),phenotype_data=phenotype_data,study_type=input$study_type,group_info_deg=paste0("group_limma",1))
      })
    
  })
  
  observeEvent(input$update_heat_kegg,{
      if(grepl("Metabolomics",input$study_type)){
        output[[paste0("heat_kegg_table",1)]] = renderDT({
          enrich_score()
        })
      }else{
        output[[paste0("heat_kegg_table",1)]] = renderDT({
          data= enrich_score()
          datatable(data, 
                    rownames = FALSE, 
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = '_all')),
                      pageLength = 10
                    )
          ) %>%
            formatSignif(columns = c('Odds.Ratio', 'Combined.Score', 'P.value', 'Adjusted.P.value'), digits = 3)
          
        })
      }
      
  })
  
  ##################################################################################################
  ###############################################Boxplot###############################################
  ##################################################################################################
  
  ########## boxplot UI#################
  output$all_boxplot_results <- renderUI({
      box_plotUI(headings=input$study_type,output_info = paste0("boxplot",1),input_info =paste0("box_gene",1),group_info = paste0("group_box",1),phenotype = phenotype_to_use())
  })
  
  observe({
      trigger = input$mainnav
      # boxplot_observe(session,input_info =paste0("box_gene",x),geneID = geneID()[[x]])
      boxplot_observe(session,input_info= paste0("box_gene",1),geneID = geneID_ID(),study_type=input$study_type)
  })


  ########## boxplot server #############
  box_res = eventReactive(input$update_box,{
      box_plotServer(input,output,session,input_gene=paste0("box_gene",1),input_group=paste0("group_box",1), data=genotype(), group_data=phenotype_to_use())
    })



  observe({
      output[[paste0("boxplot",1)]]<-renderPlotly({
        box_res()
      })
  })
  
  
  #================================================================================
  ####################### Heatmap #################################################
  #================================================================================
  
  output$all_heat_results<-renderUI({
      heatmapUI(headings = input$study_type,output_info = paste0("heat",1),input_info = paste0("gene_heat",1),group_info = paste0("group_heat",1),
                phenotype = phenotype_to_use(),color_annot=paste0("color_annot_heat",1))
      
      
  })
  
  observe({
      trigger = input$mainnav
      heatmap_gene_sel(session,input_info =paste0("gene_heat",1),genotype = genotype())
    
  })

  ###########heatmap server###########
  # heat_res = eventReactive(input$update_heat,{
  observeEvent(input$update_heat,{

      output[[paste0("heat",1)]]<-renderPlot({
        heatmapserver(input,output,session,data=genotype(),data_class=phenotype_to_use(),input_gene=paste0("gene_heat",1),group_info=paste0("group_heat",1),deg_res = deg_res(),study_type=input$study_type)
      })
    
  })
  
  #=========================================================================================
  ####################### any selected genes KEGG #################################################
  #=========================================================================================
  # Make the tab not visible for lipdomics
  tab_kegg_visible_any = reactiveVal(TRUE)
  
  observe({
    if(grepl("Lipidomics",input$study_type)|grepl("Proteomics",input$study_type)){ ################# Need to change on server
      tab_kegg_visible_any(FALSE)
    }else{
      tab_kegg_visible_any(TRUE)
    }
    
  })
  
  observe({
    if (tab_kegg_visible_any()) {
      showTab(inputId = "mainnav_kegg_any", target = "KEGG for selected genes")
    } else {
      hideTab(inputId = "mainnav_kegg_any", target = "KEGG for selected genes")
    }
  })
  
  
  output$all_keggheat_results_any<-renderUI({
    kegg_heatUI_any(headings = input$study_type,input_info =paste0("gene_heat_kegg_any",1),output_info = paste0("heat_kegg_any",1),output_info_tb = paste0("heat_kegg_table_any",1),kegg_path = paste0("kegg_path_any",1) ,group_info = paste0("group_heat_kegg_any",1),
                    phenotype = phenotype_to_use(),download_info=paste0("downloadheatkegg_any",1),gene_pathway = paste0("gene_pathway_any",1),color_annot=paste0("color_annot_any",1)
    )
    
  })
  observe({
    trigger = input$mainnav_kegg_any
    heatmap_gene_sel1(session,input_info =paste0("gene_heat_kegg_any",1),genotype = genotype(),study_type=input$study_type)
  })
  
  
  enrich_score_any<-reactive({
    enrich_res_tb_any(input,output,session,paste0("gene_heat_kegg_any",1),study_type = input$study_type)
    
  })
  
  observe({
    gene_input_id <- paste0("gene_heat_kegg_any", 1)
    observeEvent(input[[gene_input_id]], {
      # This will trigger when the gene selection changes
      enrich_res <- enrich_score_any() # fetches the current enrichment score which contains KEGG pathways
      if (is.data.frame(enrich_res)) {
        # Update the KEGG pathway dropdown
        kegg_path_choices <- enrich_res$Term
        updateSelectInput(session, paste0("kegg_path_any", 1), choices = kegg_path_choices)
      }
    }, ignoreNULL = TRUE)
    
  })
  #observeEvent(input$update_heat_kegg_any,{
  observe({
      if(grepl("Metabolomics",input$study_type)){
        output[[paste0("heat_kegg_table_any",1)]] = renderDT({
          enrich_score_any()
        })
      }else{
        output[[paste0("heat_kegg_table_any",1)]] = renderDT({
          data= enrich_score_any()
          datatable(data, 
                    rownames = FALSE, 
                    options = list(
                      columnDefs = list(list(className = 'dt-left', targets = '_all')),
                      pageLength = 10
                    )
          ) %>%
            formatSignif(columns = c('Odds.Ratio', 'Combined.Score', 'P.value', 'Adjusted.P.value'), digits = 3)
          
        })
      }
      
    
  })
  
  observe({
    
      trigger = input$mainnav_kegg_any
      enrich_res = enrich_score_any()[]
      enrich_res = data.frame(enrich_res)
      
      kegg_sel_observe(session,input_info =paste0("kegg_path_any",1),enrich_res = enrich_res)
    
  })
  

  
  #kegg_heat = eventReactive(input$update_heat_kegg,{
  observeEvent(input$update_heat_kegg_any,{
      enrich_data = enrich_score_any()
      phenotype_data = phenotype_to_use()
      omics_data = genotype()
      
      output[[paste0("heat_kegg_any",1)]] = renderPlot({
        kegg_heatServer_any(input,output,session,enrich_data =enrich_data,omics_data =omics_data ,kegg_path= paste0("kegg_path_any",1),group_info=paste0("group_heat_kegg_any",1),phenotype_data=phenotype_data,study_type=input$study_type)
      })
    
  })  
  
  #=============================================================================
  ################### Scatterplot within study ####################################
  ######ui-scatter######
  output$all_scatter_results <- renderUI({
      scatter_plotUI(headings=input$study_type,output_info = paste0("scatter",1),input_info1 =paste0("scatter_gene1_",1),input_info2 = paste0("scatter_gene2_",1),
                     group_info = paste0("group_scatter",1),download_info=paste0("downloadscatter",1),phenotype = phenotype_to_use())
    
  })
  
  observe({
      trigger = input$mainnav
      scatter_observe(session=session,input_info1 =paste0("scatter_gene1_",1),input_info2 = paste0("scatter_gene2_",1),geneID = geneID_ID(),study_type=input$study_type)
    
  })
  
  scatter_res = eventReactive(input$update_scatter,{
      scatter_plot(input,output,session,data=genotype(),group_data=phenotype_to_use(),input_group=paste0("group_scatter",1),input_gene1=paste0("scatter_gene1_",1),input_gene2=paste0("scatter_gene2_",1),study_type=input$study_type)
    
  })
  
  ######server-scatter#######
  observe({

      output[[paste0("scatter",1)]]<-renderPlotly({
        scatter_res()
      })
    
  })
  
  
  # Download figure 
  observe({
      output[[paste0("downloadscatter",1)]]<-downloadHandler(
        filename = function() {
          paste0("scatter_",input$study_type,'.png', sep = '')
          
        },
        content = function(file){
          ggsave(file, scatter_res(),type="cairo-png",width = 10, height = 10,units="cm")
        }
      )
    
  })
}

shinyApp(ui = ui, server = server)
