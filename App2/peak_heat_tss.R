peak_heat_tss_UI = function(headings,output_info){
  div(id="box1",
      column(3,
             panel(
               heading=headings, status="info",
               conditionalPanel( 
                 condition="input.update_peak_tss_heat==0",
                 br(),
                 p("Please click 'Run' button to generate figures")
               ),
               
               conditionalPanel( 
                 condition="input.update_peak_tss_heat>0",
                 br(),
                 withLoader(imageOutput(output_info),type = "html",loader = "dnaspin")
                 
               )
               
             )  
      )
      
  )
}

peak_heat_tss_server = function(input,output,session,file_folder,sample_sel){
  file = file_folder
  fig = paste0(file,"/peakTSSHeatmap/",sample_sel,".bw.png")
  
  list(src=fig,alt="fig",width=200,height=400)
}