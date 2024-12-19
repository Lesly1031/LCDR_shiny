peak_meta_UI = function(headings,output_info){
  div(id="box1",
      column(6,
             panel(
               heading=headings, status="info",
               conditionalPanel( 
                 condition="input.update_peak_meta==0",
                 br(),
                 p("Please click 'Run' button to generate figures")
               ),
               
               conditionalPanel( 
                 condition="input.update_peak_meta>0",
                 br(),
                 withLoader(imageOutput(output_info),type = "html",loader = "dnaspin")
                 
               )
               
             )  
      )
      
  )
}

peak_meta_server = function(input,output,session,file_folder,sample_sel){
  file = file_folder
  fig = paste0(file,"/metaPeak/",sample_sel,".bw.png")
  
  list(src=fig,alt="fig",width=400,height=400)
}