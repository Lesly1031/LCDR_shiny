sample_UI<-function(headings,output_info){
  
  column(6,
         panel(
           status="info",heading=paste0("Study: ",headings),
           # if(grepl("lipidomics",study_type)){
           #   tagList(
           #  # radioButtons(identify_lipi,"Collapse lipidomics class feature?",c("Yes","No"),"Yes"),
           #   withLoader( rHandsontableOutput(output_info),type="html",loader="loader5")
           #   )
           # }else{
             withLoader( rHandsontableOutput(output_info),type="html",loader="loader5")
#           }

         )
  )
}

sample_server = function(input,output,session,pheno_data){
  res = rhandsontable(pheno_data,width = 600, height = 500, stretchH = "all")
  return(res)
}