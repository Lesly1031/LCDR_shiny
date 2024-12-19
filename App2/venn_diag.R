#==========================================================================
########## Function for venn diagram ################
#==========================================================================
venn_gram = function(input,output,session,bedgr){
  unifiedlist = reduce(do.call(c,bedgr))
  fo = lapply(1:length(bedgr),function(x){
    findOverlaps(unifiedlist,bedgr[[x]])
  })
  fo = lapply(1:length(bedgr),function(x){
    queryHits(fo[[x]])
  })
  names(fo) = input$sample_sel_venn

 # mycol = ifelse(length(bedgr)==1,"#440154ff",ifelse(length(bedgr)==2,c("#440154ff",'#21908dff'),brewer.pal(length(bedgr), "Pastel2")))
  p = venn.diagram(x = fo, filename=NULL,output=TRUE)
  grid::grid.draw(p)
  
}
