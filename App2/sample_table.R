sample_UI<-function(output_info){
  
  column(6,
         panel(
           status="info",
           heading=paste0("Phenotype data"),

             withLoader( rHandsontableOutput(output_info),type="html",loader="loader5")

         )
  )
}

sample_server = function(input,output,session,pheno_data){
  res = rhandsontable(pheno_data,width = 800, height = 500, stretchH = "all")
  return(res)
}