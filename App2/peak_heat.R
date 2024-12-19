peak_heat_UI = function(headings,output_info){
  div(id="box1",
      column(3,
             panel(
             #  heading = HTML(paste0("<span style='font-size: 10px;'>", headings, "</span>")),
               heading=headings, 
               status="info",
               conditionalPanel( 
                 condition="input.update_peak_heat==0",
                 br(),
                 p("Please click 'Run' button to generate figures")
               ),
               
               conditionalPanel( 
                 condition="input.update_peak_heat>0",
                 br(),
                 withLoader(imageOutput(output_info),type = "html",loader = "dnaspin")
                 
               )

             )  
             )

      )
}

peak_heat_server = function(input,output,session,file_folder,sample_sel){
  file = file_folder
  fig = paste0(file,"/peakHeatmap/",sample_sel,".bw.png")
  list(src=fig,alt="fig",width=200,height=400)
}
