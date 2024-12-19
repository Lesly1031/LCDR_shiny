immuneDeconvUI_table<-function(headings,output_info,input_info,deconv_method,download_info){
  column(10,
         panel(
           status="info",heading=paste0("Study: ",headings),
           dropdown(
             selectInput(deconv_method,"Choose a method", choices=NULL),
             downloadButton(download_info, "Download"),
             circle = TRUE, status = "primary", icon = icon("gear"), width = "300px",
             tooltip = tooltipOptions(title = "Click to see inputs !")
           ),
           
           conditionalPanel( 
             condition="input.update_table_immuneDeconv==0",
             br(),
             p("Please click 'Run' button to generate figures")
           ),
           
           conditionalPanel( 
             condition="input.update_table_immuneDeconv>0",
             br(),
             (div(style='max-width: 100%; height: 500px; width: 90%;overflow-y: scroll;overflow-x: scroll;',
                  div(style="text-align: center;",
                      withLoader(DTOutput(output_info,width="100%",height="500px"),type="html",loader="loader5"))
                  
                  
             ))
             # DTOutput(output_info)
           )
           
           
           
         )
  )
  
}

immuneDeconvUI_method_sel_table = function(session,study_model,input_info){
  
  if(study_model=="mouse"){
    updateSelectizeInput(session,input_info,choices = c("mmcp_counter","dcq","base"),selected="base",server = T)
  }else{
    updateSelectizeInput(session,input_info,choices = c("xcell","timer","cibersort","quantiseq","cibersort_abs",
                                                        "mcp_counter","epic","abis","consensus_tme","estimate"),server = T)
  }
  
}


  