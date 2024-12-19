################# 030321 ######################
# Include considering metabolomics #
#============ 102621 ====================
# standardize RNAseq study (2808) and metabolomics (3230)
# speedup fuzzy search
# Make the checkbox only avaliable for metabolomics study 101422
################################################
packages = c("shiny","ggplot2","ggpubr","tidyverse","DT","data.table","plotly","shinyjs","shinyBS","pheatmap","DESeq2","edgeR",
             "msigdbr","shinycustomloader","jose","EBSeq","webshot","shinyWidgets","org.Mm.eg.db","grid","xCell","immunedeconv",
             "org.Hs.eg.db","ggrepel","fst","data.table","DEGreport","nipals","enrichR","stringi","rhandsontable")



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

    else {
      install.packages(x, dependencies = TRUE,repos = "http://cran.us.r-project.org")
    }
  }
  library(x, character.only = TRUE)
})


source("scatter_plot.R")
source("box_plot.R")
source("pca_plot.R")
source("DEG_table.R")
source("DEG_fig.R")
source("heat_kegg.R")
source("heatmap_plot.R")
source("heat_deg.R")
source("sample_table.R")
source("immune_deconv_table.R")
source("heat_immuneDeconv.R")
source("box_immuneDeconv.R")

source("scatter_plot_across.R")
source("heat_kegg_any.R")


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
                                                       
                                                       actionButton("update_DEG","Run"),
                                                       br(),
                                                       p('If you select an RNAseq study, you will have the option to download the normalized data following voom transformation. Simply click on the gear icon in the top left corner of the study panel and select "Download voom data." Additionally, by clicking on "Download DEG table," you can download the table of differentially expressed genes (DEG).')
                                                ),
                                                column(10,
                                                       uiOutput("all_deg_results")
                                                )
                                                
                                              )       
                                     ),
                                     tabPanel("Vocalno plot and MA plot",
                                              fluidRow(
                                                column(2,selectInput('plot_type_MA',"Select a figure type",choices = c("MA plot","Volcano plot")),
                                                       selectInput("which_p","Choose p value or adjusted p",choices = c("pval","adjust.p"),selected="pval"),
                                                       numericInput("p_threshold_MA", "pval/adj.pval threshold", 0.05,min=0.0001,max=0.1,step=0.05),
                                                       numericInput("fc_threshold_MA", "Select fold change threshold", value=1.5),
                                                       actionButton("update_MA","Run")
                                                ),
                                                column(10,uiOutput("all_deg_figs"))
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
                
                tabPanel("Heatmap",
                         tabsetPanel(type="tabs",id="mainnav_kegg_any",
                           tabPanel("Heatmap",
                                    fluidRow(
                                      column(2,
       
                                           checkboxInput("show_row_names",
                                                           "Show omics names", value = FALSE),

                                             p("Shows the top 25 markers from DEG results by default"),
                                             actionButton("update_heat","Run")
                                             
                                      ),
                                      column(10,uiOutput("all_heat_results"))
                                    )
                           ),
                           tabPanel("KEGG for selected genes",
                                    fluidRow(
                                      column(2,

                                             checkboxInput("show_row_names_kegg_any",
                                                           "Show gene names", value = FALSE),
                                           
                                             p("Click the button below after selecting the pathway from the drop down gear button"),
                                             actionButton("update_heat_kegg_any","Run")
                                      ),
                                      column(10,uiOutput("all_keggheat_results_any"))
                                    )
                           )
             
                           
                           
                         )
                         
                         
                         
                         
                ),
                tabPanel("Scatterplot",
                         tabsetPanel(type="tabs",
                                     tabPanel("For each study",
                                              fluidRow(
                                                column(2,
                                                       actionButton("update_scatter","Run")
                                                       
                                                ),
                                                column(10,
                                                       uiOutput("all_scatter_results"),
                                                       textOutput("text")
                                                )
                                                
                                              )
                                     )
                                    
                         )
                         
                ),
                
                
                
                
                
                

                tabPanel("Immune decovolution",
                         tabsetPanel(type="tabs",id="mainnav_immunedeconv", 
                           tabPanel("Immune deconv table",
                                    fluidRow(
                                      column(2,
                                             p("This tab displays heatmaps for immune deconvolution in RNA-seq studies. The data has been transformed using z-scores"),
                                             br(),
                                             p("Click 'Run' button to get the table after selection"),
                                             actionButton("update_table_immuneDeconv","Run"),
                                             
                                      ),
                                      column(10,
                                             uiOutput("all_immuneDeconv_table"),
                                             DTOutput("test_deconv")
                                      )
                                    )
                                    ),
                           tabPanel("Immune deconv heatmap",
                                    fluidRow(
                                      column(2,
                                             br(),
                                             p("Click 'Run' button to get the heatmap after selection"),
                                             actionButton("update_heat_immuneDeconv","Run"),
                                             
                                      ),
                                      column(10,
                                             uiOutput("all_immuneDeconv_heat")
                                      )
                                    )
                                    ),
                           tabPanel("immune deconv boxplot",
                                    fluidRow(
                                      column(2,
                                             br(),
                                             p("Click 'Run' button to get the boxplot after selection"),
                                             actionButton("update_box_immuneDeconv","Run"),
                                             
                                      ),
                                      column(10,
                                             uiOutput("all_immuneDeconv_box")
                                      )
                                    )
                                    )
                         )

                         )
                
    )
  )
  
  
  
)


server <- function(input, output,session) { 

  
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
  
  
  folder_direct <- reactive({
    file =  sapply(1:length(study_id()),function(x){
      lapply(1:length(study_type()[[x]]),function(y){
        file =  paste0("studies/",study_id()[x])
        paste0(file,"/",study_type()[[x]][y])
      })
    })
    unlist(file)
  })
  
  
  
  headings = reactive({
    file =  sapply(1:length(study_id()),function(x){
      lapply(1:length(study_type()[[x]]),function(y){
        paste0("Study",study_id()[[x]],"_",study_type()[[x]][y])
      })
    })
    unlist(file)
  })
  
  #===============================================================================
  ########################## phenotype data file ##############################
  #===============================================================================
  output$all_sample<-renderUI({
    sample_list <- lapply(1:length(folder_direct()),function(x){
      sample_UI(headings = headings()[[x]],output_info = paste0("sample",x))
    })
    tagList(sample_list)
  })
  
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("sample",x)]]<-renderRHandsontable({
        sample_server(input,output,session,pheno_data = phenotype()[[x]])
      })
    })
  })
  
  phenotype_to_use <- reactive({
    lapply(1:length(folder_direct()),function(x){
      dt_class=hot_to_r(input[[paste0("sample",x)]])
      dt_class=data.frame(dt_class,check.names = F)
      dt_class = dt_class[complete.cases(dt_class$X),]
      dt_class
    })
  })
  
  genotype <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics.csv")),check.names = F)
    })
  })
  genotype_raw <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics_raw.csv")),check.names = F)
    })
  })
  
  genotype_identify <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics_identify.csv")),check.names = F)
    })
  })
  genotype_identify_raw <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics_identify_raw.csv")),check.names = F)
    })
  })
  
  # for lipidomics class colllapse
  genotype_class <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics_class.csv")),check.names = F)
    })
  })
  genotype_raw_class <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics_raw_class.csv")),check.names = F)
    })
  })

  genotype_identify_class <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics_identify_class.csv")),check.names = F)
    })
  })
  genotype_identify_raw_class <- reactive({
    lapply(1:length(folder_direct()),function(x){
      data.frame(fread(paste0(folder_direct()[x],"/omics_identify_raw_class.csv")),check.names = F)
    })
  })
  
  geneID_class <- reactive({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav_sub
      genotype_class()[[x]][,"ID"]

    })
  })
  
  
  
  gene_matrix_immuneDeconv <- reactive({
    lapply(1:length(folder_direct()),function(x){
      dt = data.frame(fread(paste0(folder_direct()[x],"/gene_matrix_deconv.csv")))
      rownames(dt) = dt$V1
      dt = data.frame(dt[,-1])
      dt
    })
  })
  
  
  
  phenotype <- reactive({
    lapply(1:length(folder_direct()),function(x){
      dt = data.frame(fread(paste0(folder_direct()[x],"/phenotype_data.csv")),check.names = F)
      dt[dt==""] = NA
      dt
    })
  })
  
  immune_decov_matric <- reactive({
    tryCatch({
      lapply(1:length(folder_direct()),function(x){
        dt = fread(paste0(folder_direct()[x],"/",input[[paste0("deconv_method_table",x)]],"_deconv.csv"))
        
      })
    },error = function(e){
      ""
    })

  })
  

  
  geneID <- reactive({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav_sub
      #      genotype()[[x]][,"ID"]
      # if(input$pca_meta=="Yes"|input$deg_meta=="Yes"){
      #   genotype_identify()[[x]][,"ID"]
      # }else{
      genotype()[[x]][,"ID"]
      #rownames(deg_res()[[x]])
      #}
    })
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
  
  #===============================================================================
  ######################### phenotype data file #############################
  #===============================================================================
  output$all_sample<-renderUI({
    sample_list <- lapply(1:length(folder_direct()),function(x){
      sample_UI(headings = headings()[[x]],output_info = paste0("sample",x))
    })
    tagList(sample_list)
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("sample",x)]]<-renderRHandsontable({
        sample_server(input,output,session,pheno_data = phenotype()[[x]])
      })
    })
  })
  
  phenotype_to_use <- reactive({
    lapply(1:length(folder_direct()),function(x){
      dt_class=hot_to_r(input[[paste0("sample",x)]])
      dt_class=data.frame(dt_class,check.names = F)
      dt_class = dt_class[complete.cases(dt_class$X),]
      dt_class
    })
  })
  
  # output$test_sample = renderDT({
  #   phenotype_to_use()[[1]]
  # }
  # )
  
  #==============================================================================
  ############################# PCA ############################################
  #==============================================================================
  
  ############################PCA#############################
  output$all_pca_results<-renderUI({
    pca_list <- lapply(1:length(folder_direct()),function(x){
      # genotype = values[[paste0("genotype",x)]]
      pca_UI(headings = headings()[[x]],output_info = paste0("pca",x),group_info1 = paste0("group_pca1",x),group_info2 = paste0("group_pca2",x),download_info=paste0("downloadpca",x),
             phenotype = phenotype_to_use()[[x]],identify_meta = paste0("pca_meta",x),study_type=study_type()[[x]])
    })
    tagList(pca_list)
  })
  
  ###########pca server###########
  pca_res = eventReactive(input$update_pca,{
    lapply(1:length(folder_direct()),function(x){
      pca_server(input,output,session,data_identify=genotype_identify()[[x]],data=genotype()[[x]],data_class=phenotype_to_use()[[x]],data_rna = genotype_raw()[[x]],
                 group_info1=paste0("group_pca1",x),group_info2=paste0("group_pca2",x),folder_direct=folder_direct()[[x]],identify_meta = paste0("pca_meta",x))
    })
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("pca",x)]]<-renderPlotly({
        p =  pca_res()[[x]]
        ggplotly(p, tooltip="text") %>%
          layout(legend = list(orientation = 'h', x = 1.02, y = 1))
      })
    })
  })
  
  # Download figure 
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("downloadpca",x)]]<-downloadHandler(
        filename = function() {
          paste0("PCA_",headings()[[x]],'.png', sep = '')
          
        },
        content = function(file){
          ggsave(file, pca_res()[[x]],type="cairo-png",width = 15, height = 15,units="cm")
        }
      )
    })
  })
  
  ########################### DEG table#################################################
  
  
  
  output$all_deg_results<-renderUI({
    deg_list <- lapply(1:length(folder_direct()),function(x){
      # genotype = values[[paste0("genotype",x)]]
      limma_UI(headings=headings()[[x]],output_info = paste0("limma",x),voom_dt = paste0("voom",x),group_info = paste0("group_limma",x),
               download_info=paste0("downloadlimma",x),download_info2=paste0("downloadvoom",x),
               phenotype = phenotype_to_use()[[x]], 
               deg_meta = paste0("deg_meta",x),
               study_type = study_type()[[x]])
    })
    tagList(deg_list)
  })
  
  
  
  
  ######## DEG table server ##########
  deg_res <- eventReactive(input$update_DEG,{
    lapply(1:length(folder_direct()),function(x){
      limma_server(input,output,session,folder_direct=folder_direct()[[x]],data_identify = genotype_identify_raw()[[x]],data_lip=genotype_identify_raw_class()[[x]],data=genotype_raw()[[x]],data_class = phenotype_to_use()[[x]],
                   group_info = paste0("group_limma",x),study_type =study_type()[[x]],deg_meta = paste0("deg_meta",x))
    })
  })
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("limma",x)]]<-renderDT({
        deg_res = deg_res()[[x]]
        # deg_res$rowID =  sapply( str_split(rownames(deg_res), "_"), "[", 1)
        deg_res$rowID = rownames(deg_res)
        rownames(deg_res) = NULL
        
        #deg_res = data.frame(deg_res[,c("rowID",names(deg_res)[-ncol(deg_res)]  )])
        
        #        if(grepl("metabolomics",folder_direct()[[x]])|grepl("lipidomics",folder_direct()[[x]])){
        if(grepl("metabolomics",folder_direct()[[x]])){
          
          # names = sapply(strsplit(deg_res$rowID,"_"),"[",1) 
          # IDs = unlist(lapply(1:length(deg_res$rowID),function(x){
          #   paste(strsplit(deg_res$rowID[x],"_")[[1]][2:3],collapse = "_")
          # }))
          # deg_res$Names= ifelse(names=="X",deg_res$rowID,names)
          # deg_res$IDs= ifelse(names=="X",deg_res$rowID,IDs)
          split_sample_names <- function(names) {
            # The first part up to the polarity (pos or neg)
            part1 <- str_extract(names, ".*(?=_pos|_neg)")
            # The second part starting from the polarity (pos or neg)
            part2 <- str_extract(names, "(pos|neg)_.+")
            # Combine the parts into a data frame
            data.frame(part1, part2)
          }
          
          # Apply the function to split the names
          split_names_df <- split_sample_names(deg_res$rowID)
          names = split_names_df$part1
          IDs = split_names_df$part2
          deg_res$Names= ifelse(names=="X",deg_res$rowID,names)
          deg_res$IDs= ifelse(names=="X",deg_res$rowID,IDs)
          
        }else if(grepl("lipidomics",folder_direct()[[x]])){
          deg_res$Names = deg_res$rowID
          deg_res$IDs = ""
        }else{
          deg_res$Names = sapply(strsplit(deg_res$rowID,"_"),"[",1)
          deg_res$IDs = sapply(strsplit(deg_res$rowID,"_"),"[",-1)
          
        }
        
        deg_res = data.frame(deg_res[,c("Names","IDs","log2FoldChange","baseMean","pvalue","padj")],check.names=F)
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
  })
  
  ########################### download DEG tables ################################
  # Download TABLE 
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("downloadlimma",x)]]<-downloadHandler(
        filename = function() {
          paste0("DEG_",headings()[[x]],'.csv', sep = '')
          
        },
        content = function(file){
          fwrite(deg_res()[[x]],file,row.names = T)
        }
      )
    })
  })
  
  ######################## DEG voom results #################################
  voom_res <- eventReactive(input$update_DEG,{
    
    lapply(1:length(folder_direct()),function(x){
      if(grepl("rna-seq",folder_direct()[[x]])){ 
        limma_voom(input,output,session,folder_direct=folder_direct()[[x]],data=genotype_raw()[[x]],data_class = phenotype_to_use()[[x]],
                   group_info = paste0("group_limma",x))
      }
    })
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      if(grepl("rna-seq",folder_direct()[[x]])){ 
        
        output[[paste0("voom",x)]]<-renderDT({
          res = voom_res()[[x]]
          datatable(
            res ,
            rownames = T,plugins = "ellipsis",
            extensions = 'Buttons',
            options = list(
              autoWidth = FALSE, scrollX = TRUE,
              columnDefs = list(list(
                width = "10px", targets = "_all",
                render = JS("$.fn.dataTable.render.ellipsis( 30, false )")
              )))
          )
          
        })
      }
    })
  })
  
  # Download TABLE 
  observe({
    lapply(1:length(folder_direct()),function(x){
      if(grepl("rna-seq",folder_direct()[[x]])){ 
        
        output[[paste0("downloadvoom",x)]]<-downloadHandler(
          filename = function() {
            paste0("voom_",headings()[[x]],'.csv', sep = '')
            
          },
          content = function(file){
            fwrite(voom_res()[[x]],file,row.names = T)
          }
        )
        
      }
    })
  })
  
  ##############################DEG visualization ################################
  
  ########## deg figure server #############
  output$all_deg_figs<-renderUI({
    deg_fig_list <- lapply(1:length(folder_direct()),function(x){
      # genotype = values[[paste0("genotype",x)]]
      deg_figUI(headings = headings()[[x]],output_info = paste0("deg_fig",x),
                download_info=paste0("downloaddegfig",x),fill=paste0("deg_gene",x))
    })
    tagList(deg_fig_list)
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav_sub
      deg_fig_observe(session=session,input_info =paste0("deg_gene",x),deg_res = deg_res()[[x]])
    })
  })
  ########## DEG visualize server ############
  
  MA_limma<-eventReactive(input$update_MA,{
    lapply(1:length(folder_direct()),function(x){
      deg_figServer(input,output,session,data.frame(deg_res()[[x]]),folder_direct=folder_direct()[[x]],fill_ma =paste0("deg_gene",x),title = input[[paste0("group_limma",x)]] )
    })
    
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("deg_fig",x)]]<-renderPlotly({
        ggplotly(MA_limma()[[x]],tooltip="text")
        
      })
    })
  })
  
  # Download figure 
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("downloaddegfig",x)]]<-downloadHandler(
        filename = function() {
          paste0("DEG_fig_",headings()[[x]],'.png', sep = '')
          
        },
        content = function(file){
          ggsave(file, MA_limma()[[x]],type="cairo-png",width = 25, height = 15,units="cm")
        }
      )
    })
  })
  
  
  ####################### Heatmap DEG ######################################
  
  output$all_heat_deg<-renderUI({
    heat_list <- lapply(1:length(folder_direct()),function(x){
      heatmapDegUI(headings = headings()[[x]],output_info = paste0("heat_deg",x),color_annot=paste0("color_annot_heat_deg",x))
      
      
    })
    tagList(heat_list)
  })
  
  observeEvent(input$update_heat_deg,{
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("color_annot_heat_deg",x)]]<-renderText({
        heatDeg_color_annot(input,output,session, study_type=study_type()[[x]])
      })
    })
  })
  
  
  
  ###########heatmap server###########
  # heat_res = eventReactive(input$update_heat,{
  observeEvent(input$update_heat_deg,{
    lapply(1:length(folder_direct()),function(x){
      
      output[[paste0("heat_deg",x)]]<-renderPlot({
        heatmapDeg_server(input,output,session,data=genotype()[[x]],data_lip=genotype_class()[[x]],data_identify=genotype_identify()[[x]],data_class=phenotype_to_use()[[x]],group_info=paste0("group_limma",x),deg_res = deg_res()[[x]],folder_direct = folder_direct()[[x]],deg_meta= paste0("deg_meta",x),group_info_deg=paste0("group_limma",x))
      })
    })
  })
  
  
  
  # output$test_degheat = renderDT({
  #  
  #   heatmapDeg_server(input,output,session,data=genotype()[[1]],data_class=phenotype_to_use()[[1]],
  #                     group_info=paste0("group_limma","1"),deg_res = deg_res()[[1]])
  #   
  # })
  
  ###########################################################################################
  #==============================================================================================
  ##################################################################################################
  ###############################################Boxplot###############################################
  ##################################################################################################
  
  ########## boxplot UI#################
  output$all_boxplot_results <- renderUI({
    box_list <- lapply(1:length(folder_direct()),function(x){
      box_plotUI(headings=headings()[[x]],output_info = paste0("box",x),input_info =paste0("box_gene",x),group_info = paste0("group_box",x),phenotype = phenotype_to_use()[[x]])
    })
    tagList(box_list)
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav
      
      # boxplot_observe(session,input_info =paste0("box_gene",x),geneID = geneID()[[x]])
      boxplot_observe(session,input_info= paste0("box_gene",x), geneID()[[x]],geneID_lip=geneID_class()[[x]],study_type=study_type()[[x]])
    })
  })
  
  ########## boxplot server #############
  box_res = eventReactive(input$update_box,{
    lapply(1:length(folder_direct()),function(x){
      box_plotServer(input,output,session,input_gene=paste0("box_gene",x),input_group=paste0("group_box",x), data=genotype()[[x]], group_data=phenotype_to_use()[[x]],folder_direct=folder_direct()[[x]],study_type=study_type()[[x]],data_lip=genotype_class()[[x]],data_rnaseq = genotype_raw()[[x]])
    })
  })
  
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("box",x)]]<-renderPlotly({
        box_res()[[x]]
      })
    })
  })
  
  
  
  #================================================================================
  ####################### Heatmap #################################################
  #================================================================================
  
  output$all_heat_results<-renderUI({
    heat_list <- lapply(1:length(folder_direct()),function(x){
      heatmapUI(headings = headings()[[x]],output_info = paste0("heat",x),input_info = paste0("gene_heat",x),group_info = paste0("group_heat",x),
                phenotype = phenotype_to_use()[[x]],color_annot=paste0("color_annot_heat",x))
      
      
    })
    tagList(heat_list)
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav
      heatmap_gene_sel(session,input_info =paste0("gene_heat",x),genotype = genotype()[[x]],genotype_class = genotype_class()[[x]],study_type=study_type()[[x]])
    })
  })
  
  observeEvent(input$update_heat,{
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("color_annot_heat",x)]]<-renderText({
        heat_color_annot(input,output,session, study_type=study_type()[[x]])
      })
    })
  })
  
  
  
  ###########heatmap server###########
  # heat_res = eventReactive(input$update_heat,{
  observeEvent(input$update_heat,{
    lapply(1:length(folder_direct()),function(x){
      
      output[[paste0("heat",x)]]<-renderPlot({
        heatmapserver(input,output,session,data=genotype()[[x]],data_class=phenotype_to_use()[[x]],input_gene=paste0("gene_heat",x),group_info=paste0("group_heat",x),deg_res = deg_res()[[x]],folder_direct = folder_direct()[[x]],data_lip=genotype_class()[[x]],study_type=study_type()[[x]])
      })
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
    lapply(1:length(folder_direct()),function(x){
      if(grepl("lipidomics",folder_direct()[[x]])|grepl("proteomics",folder_direct()[[x]])){ ################# Need to change on server
        tab_kegg_visible(FALSE)
      }else{
        tab_kegg_visible(TRUE)
      }
    })
    
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
    heat_list_kegg <- lapply(1:length(folder_direct()),function(x){
      kegg_heatUI(headings = headings()[[x]],output_info = paste0("heat_kegg",x),output_info_tb = paste0("heat_kegg_table",x),kegg_path = paste0("kegg_path",x) ,group_info = paste0("group_heat_kegg",x),
                  phenotype = phenotype_to_use()[[x]],download_info=paste0("downloadheatkegg",x),gene_pathway = paste0("gene_pathway",x),color_annot=paste0("color_annot",x)
      )
    })
    tagList(heat_list_kegg)
  })
  
  enrich_score<-reactive({
    lapply(1:length(folder_direct()),function(x){
      enrich_res_tb(input,output,session,deg_data = data.frame(deg_res()[[x]]),folder_direct =folder_direct()[[x]])
    })
    
  })
  
  
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav_sub
      enrich_res = enrich_score()[[x]]
      #    if(input$omics_filter_kegg=="No"){
      enrich_res = data.frame(enrich_res)
      # }else{
      #   enrich_res = subset(enrich_res,enrich_res$length_genes<=input$number_omics_kegg&enrich_res$length_genes>=3)
      # }
      
      kegg_sel_observe(session,input_info =paste0("kegg_path",x),enrich_res = enrich_res)
    })
  })
  
  # observeEvent(input$update_heat_kegg,{
  observe({  
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("color_annot",x)]]<-renderText({
        kegg_color_annot(input,output,session, folder_direct=folder_direct()[[x]])
      })
    })
  })
  
  #kegg_heat = eventReactive(input$update_heat_kegg,{
  observeEvent(input$update_heat_kegg,{
    
    lapply(1:length(folder_direct()),function(x){
      enrich_data = enrich_score()[[x]]
      phenotype_data = phenotype_to_use()[[x]]
      omics_data = genotype_identify_raw()[[x]]
      
      output[[paste0("heat_kegg",x)]] = renderPlot({
        kegg_heatServer(input,output,session,enrich_data =enrich_data,omics_data =omics_data ,kegg_path= paste0("kegg_path",x),group_info=paste0("group_heat_kegg",x),phenotype_data=phenotype_data,folder_direct = folder_direct()[[x]],group_info_deg=paste0("group_limma",x))
      })
    })
  })
  
  observeEvent(input$update_heat_kegg,{
    lapply(1:length(folder_direct()),function(x){
      if(grepl("metabolomics",folder_direct()[[x]])){
        output[[paste0("heat_kegg_table",x)]] = renderDT({
          enrich_score()[[x]]
        })
      }else{
        output[[paste0("heat_kegg_table",x)]] = renderDT({
          data= enrich_score()[[x]]
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
  })
  
  #===============================================================================
  ###################### Scatterplot ###########################
  #=============================================================================
  ################### Scatterplot within study ####################################
  ######ui-scatter######
  output$all_scatter_results <- renderUI({
    scatter_list <- lapply(1:length(folder_direct()),function(x){
      scatter_plotUI(headings=headings()[[x]],output_info = paste0("scatter",x),input_info1 =paste0("scatter_gene1_",x),input_info2 = paste0("scatter_gene2_",x),
                     group_info = paste0("group_scatter",x),download_info=paste0("downloadscatter",x),phenotype = phenotype_to_use()[[x]])
    })
    tagList(scatter_list)
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav
      scatter_observe(session=session,input_info1 =paste0("scatter_gene1_",x),input_info2 = paste0("scatter_gene2_",x),geneID = geneID()[[x]],geneID_lip=geneID_class()[[x]],study_type=study_type()[[x]])
    })
  })
  
  scatter_res = eventReactive(input$update_scatter,{
    lapply(1:length(folder_direct()),function(x){ 
      scatter_plot(input,output,session,data=genotype()[[x]],group_data=phenotype_to_use()[[x]],input_group=paste0("group_scatter",x),input_gene1=paste0("scatter_gene1_",x),input_gene2=paste0("scatter_gene2_",x),study_type=study_type()[[x]],data_lip=genotype_class()[[x]])
    })
  })
  
  ######server-scatter#######
  observe({
    lapply(1:length(folder_direct()),function(x){
      
      output[[paste0("scatter",x)]]<-renderPlotly({
        scatter_res()[[x]]
      })
    })
  })
  
  
  # Download figure 
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("downloadscatter",x)]]<-downloadHandler(
        filename = function() {
          paste0("scatter_",headings()[[x]],'.png', sep = '')
          
        },
        content = function(file){
          ggsave(file, scatter_res()[[x]],type="cairo-png",width = 10, height = 10,units="cm")
        }
      )
    })
  })
  
  
  #--------------------------------------------------------------------------------
  ##################### Immune decovolution for RNAseq study ####################
  #---------------------------------------------------------------------------------
  # Make the tab visible for RNAseq studies only
  tab_immune_visible = reactiveVal(TRUE)
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      if(grepl("rna-seq",folder_direct()[[x]])){ ################# Need to change on server
        tab_immune_visible(TRUE)
      }else{
        tab_immune_visible(FALSE)
      }
    })
    
  })
  
  observe({
    if (tab_immune_visible()) {
      showTab(inputId = "mainnav", target = "Immune decovolution")
    } else {
      hideTab(inputId = "mainnav", target = "Immune decovolution")
    }
  })
  #---------------------------------------------------------------------------
  #################### Immune deconvolution table ###########################
  #---------------------------------------------------------------------------
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav
      immuneDeconvUI_method_sel_table(session,study_model = study_model()[[x]],input_info =paste0("deconv_method_table",x))
    })
  })
  
  output$all_immuneDeconv_table<-renderUI({
    immuneDeconv_list <- lapply(1:length(folder_direct()),function(x){
      immuneDeconvUI_table(headings = headings()[[x]],output_info = paste0("immuneDeconv_table_",x),deconv_method = paste0("deconv_method_table",x),
                           download_info = paste0("deconv_table_download",x)
      )
    })
    tagList(immuneDeconv_list)
  })
  
  ######## immune deconv table server ##########
  
  
  immune_deconv_table <- eventReactive(input$update_table_immuneDeconv,{
    lapply(1:length(folder_direct()),function(x){
      immune_decov_matric()[[x]]
    })
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[ paste0("immuneDeconv_table_",x)]]<-renderDT({
        immune_deconv_table()[[x]]
      })
    })
  })
  
  ########################### download immune deconv tables ################################
  # Download TABLE 
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("deconv_table_download",x)]]<-downloadHandler(
        filename = function() {
          paste0("immune_deconv_table_",headings()[[x]],'.csv', sep = '')
          
        },
        content = function(file){
          fwrite(immune_deconv_table()[[x]],file,row.names = T)
        }
      )
    })
  })
  
  #---------------------------------------------------------------------------
  #################### Immune deconvolution heatmap ###########################
  #---------------------------------------------------------------------------
  output$all_immuneDeconv_heat<-renderUI({
    immuneDeconv_list <- lapply(1:length(folder_direct()),function(x){
      immuneDeconvUI_heat(headings = headings()[[x]],output_info = paste0("immuneDeconv_heat_",x),group_info = paste0("group_immuneDeconv_heat",x),phenotype = phenotype_to_use()[[x]]
                          #   download_info = paste0("deconv_table_download",x)
      )
    })
    tagList(immuneDeconv_list)
  })
  
  observeEvent(input$update_heat_immuneDeconv,{
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("immuneDeconv_heat_",x)]]<-renderPlot({
        immuneDeconv_heat_server(input,output,session,data=immune_deconv_table()[[x]],data_class=phenotype_to_use()[[x]],group_info=paste0("group_immuneDeconv_heat",x))
      })
    })
  })
  
  #-----------------------------------------------------------------------------
  ###################### Immune deconv boxplot ###############################
  #-----------------------------------------------------------------------------
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav_immunedeconv
      immuneDeconvUI_method_sel_box(session,data=immune_deconv_table()[[x]],cell_type =paste0("cell_type_box",x))
    })
  })
  
  output$all_immuneDeconv_box<-renderUI({
    immuneDeconv_list <- lapply(1:length(folder_direct()),function(x){
      immuneDeconvUI_box(headings = headings()[[x]],output_info = paste0("immuneDeconv_box_",x),
                         group_info = paste0("group_immuneDeconv_box",x),phenotype = phenotype_to_use()[[x]], cell_type = paste0("cell_type_box",x)
      )
    })
    tagList(immuneDeconv_list)
  })
  
  ########## immune Deconv boxplot server #############
  immuneDconv_box_res = eventReactive(input$update_box_immuneDeconv,{
    lapply(1:length(folder_direct()),function(x){
      immuneDeconv_box_server(input,output,session,input_gene=paste0("cell_type_box",x),input_group=paste0("group_immuneDeconv_box",x),
                              data=immune_deconv_table()[[x]],group_data=phenotype_to_use()[[x]],folder_direct=folder_direct()[[x]])
    })
  })
  
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("immuneDeconv_box_",x)]]<-renderPlotly({
        immuneDconv_box_res()[[x]]
      })
    })
  })
  
  
  #=========================================================================================
  ####################### any selected genes KEGG #################################################
  #=========================================================================================
  # Make the tab not visible for lipdomics
  tab_kegg_visible_any = reactiveVal(TRUE)
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      if(grepl("lipidomics",folder_direct()[[x]])|grepl("proteomics",folder_direct()[[x]])){ ################# Need to change on server
        tab_kegg_visible_any(FALSE)
      }else{
        tab_kegg_visible_any(TRUE)
      }
    })
    
  })
  
  observe({
    if (tab_kegg_visible_any()) {
      showTab(inputId = "mainnav_kegg_any", target = "KEGG for selected genes")
    } else {
      hideTab(inputId = "mainnav_kegg_any", target = "KEGG for selected genes")
    }
  })
  
  
  output$all_keggheat_results_any<-renderUI({
    heat_list_kegg <- lapply(1:length(folder_direct()),function(x){
      kegg_heatUI_any(headings = headings()[[x]],input_info =paste0("gene_heat_kegg_any",x),output_info = paste0("heat_kegg_any",x),output_info_tb = paste0("heat_kegg_table_any",x),kegg_path = paste0("kegg_path_any",x) ,group_info = paste0("group_heat_kegg_any",x),
                      phenotype = phenotype_to_use()[[x]],download_info=paste0("downloadheatkegg_any",x),gene_pathway = paste0("gene_pathway_any",x),color_annot=paste0("color_annot_any",x)
      )
    })
    tagList(heat_list_kegg)
  })
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav_kegg_any
      heatmap_gene_sel1(session,input_info =paste0("gene_heat_kegg_any",x),genotype = genotype()[[x]],genotype_class = genotype_class()[[x]],study_type=study_type()[[x]])
    })
  })
  
  enrich_score_any<-reactive({
    lapply(1:length(folder_direct()),function(x){
      enrich_res_tb_any(input,output,session,paste0("gene_heat_kegg_any",x),folder_direct =folder_direct()[[x]])
    })
    
  })
  
  observe({
    lapply(1:length(folder_direct()), function(x) {
      gene_input_id <- paste0("gene_heat_kegg_any", x)
      observeEvent(input[[gene_input_id]], {
        # This will trigger when the gene selection changes
        enrich_res <- enrich_score_any()[[x]] # fetches the current enrichment score which contains KEGG pathways
        if (is.data.frame(enrich_res)) {
          # Update the KEGG pathway dropdown
          kegg_path_choices <- enrich_res$Term
          updateSelectInput(session, paste0("kegg_path_any", x), choices = kegg_path_choices)
        }
      }, ignoreNULL = TRUE)
    })
  })
  
  #observeEvent(input$update_heat_kegg_any,{
  observe({
    lapply(1:length(folder_direct()),function(x){
      if(grepl("metabolomics",folder_direct()[[x]])){
        output[[paste0("heat_kegg_table_any",x)]] = renderDT({
          enrich_score_any()[[x]]
        })
      }else{
        output[[paste0("heat_kegg_table_any",x)]] = renderDT({
          data= enrich_score_any()[[x]]
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
  })
  
  observe({
    lapply(1:length(folder_direct()),function(x){
      trigger = input$mainnav_kegg_any
      enrich_res = enrich_score_any()[[x]]
      enrich_res = data.frame(enrich_res)
      
      kegg_sel_observe(session,input_info =paste0("kegg_path_any",x),enrich_res = enrich_res)
    })
  })
  
  observeEvent(input$update_heat_kegg_any,{
    lapply(1:length(folder_direct()),function(x){
      output[[paste0("color_annot_any",x)]]<-renderText({
        kegg_color_annot(input,output,session, folder_direct=folder_direct()[[x]])
      })
    })
  })
  
  #kegg_heat = eventReactive(input$update_heat_kegg,{
  observeEvent(input$update_heat_kegg_any,{
    lapply(1:length(folder_direct()),function(x){
      enrich_data = enrich_score_any()[[x]]
      phenotype_data = phenotype_to_use()[[x]]
      omics_data = genotype_identify_raw()[[x]]
      
      output[[paste0("heat_kegg_any",x)]] = renderPlot({
        kegg_heatServer_any(input,output,session,enrich_data =enrich_data,omics_data =omics_data ,kegg_path= paste0("kegg_path_any",x),group_info=paste0("group_heat_kegg_any",x),phenotype_data=phenotype_data,folder_direct = folder_direct()[[x]])
      })
    })
  })
}

shinyApp(ui = ui, server = server)
