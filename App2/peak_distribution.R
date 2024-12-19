peak_distribution_UI = function(headings,output_info1,output_info2,output_download1,output_download2){
  column(11,
         panel(
           heading=headings, 
           status="info",
           fluidRow(
             column(6,
                    conditionalPanel( 
                      condition="input.update_peak_distribution==0",
                      br(),
                      p("Please click 'Run' button to generate figures")
                    ),
                    
                    conditionalPanel( 
                      condition="input.update_peak_distribution>0",
                      br(),
                      downloadButton(output_download1,"Download the peak distribution data"),
                      withLoader(plotOutput(output_info1),type = "html",loader = "dnaspin")
                      
                    )
             ),
             column(6,
                    conditionalPanel( 
                      condition="input.update_peak_distribution==0",
                      br(),
                      p("")
                    ),
                    conditionalPanel( 
                      condition="input.update_peak_distribution>0",
                      br(),
                      downloadButton(output_download2,"Download TF-binding loci data near TSS"),
                      plotOutput(output_info2)
                    )
                    
              
             )
           )
         )
   
  )
}

peakAnnot_data1 = function(input,output,session,file_folder,input_sample){
  load(paste0(file_folder,"peakDistribution/",input_sample,".narrowPeak.RData"))
   peakAnno@annoStat
}

peakAnnot_data2 = function(input,output,session,file_folder,input_sample){
  load(paste0(file_folder,"peakDistribution/",input_sample,".narrowPeak.RData"))
  p = plotDistToTSS(peakAnno)
  p$data
}


peakAnnot_fig1= function(input,output,session,file_folder,input_sample){
  load(paste0(file_folder,"peakDistribution/",input_sample,".narrowPeak.RData"))
   plotAnnoBar(peakAnno)
}

peakAnnot_fig2= function(input,output,session,file_folder,input_sample){
  load(paste0(file_folder,"peakDistribution/",input_sample,".narrowPeak.RData"))
  plotDistToTSS(peakAnno)
}


