#' Lineage Figures - Lineage Tree
#'
#' This group of Lineage Figure functions produce figures to assist in interpretation of lineage data.
#'
#' This function produces a lineage tree with a lineage bar to visualise the phylogenetic positions of the lineages.
#' The outputs of the seq_designation, node_info and lineage_info functions are required, along with the
#' phylogenetic tree and corresponding metadata file used as input for the sequence data, node info and lineage
#' information.
#'
#' @param lineage_info The output of the lineage_info function
#' @param node_data The output of the node_info function
#' @param tree A phylogenetic tree
#' @param metadata The metadata corresponding to the sequences in the tree, including "ID" "assignment" "country" and "year"
#' @param sequence_data The output of the seq_designation function
#' @return A lineage tree with a lineage bar to visualise the phylogenetic positions of the lineages
#' @export
lineage_tree <- function(lineage_info, node_data, tree, metadata, sequence_data) {
  tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
  tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)

  lineage_info$colour<-NA

  Colours<-c("Reds","Purples","YlOrBr","PuBuGn","YlOrRd","OrRd","PuBu","Pastel1","Greens","Greys",
             "GnBu","BuGn","RdPu","Oranges","BuPu","YlGn","PuRd","YlGnBu")

  lineages<-data.frame(lineage = lineage_info$lineage, subclade = NA)

  for (i in 1:length(lineages$lineage)) {
    lineages$subclade[i]<-strsplit(lineages$lineage[i], "_")[[1]][1]
  }

  letters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N",
               "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z","AA","AB","AC","AD","AE","AF","AG"
               ,"AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU","AV","AW","AX","AY","AZ",
               "BA","BB","BC","BD","BE","BF","BG","BH","BI","BJ","BK","BL","BM","BN","BO","BP","BQ","BR","BS","BT"
               ,"BU","BV", "BW", "BX", "BY", "BZ","CA", "CB", "CC", "CD", "CE", "CF", "CG", "CH", "CI", "CJ", "CK", "CL", "CM", "CN",
               "CO", "CP", "CQ", "CR", "CS", "CT", "CU", "CV", "CW", "CX", "CY", "CZ","DA", "DB", "DC", "DD", "DE", "DF", "DG", "DH", "DI", "DJ", "DK", "DL", "DM", "DN",
               "DO", "DP", "DQ", "DR", "DS", "DT", "DU", "DV", "DW", "DX", "DY", "DZ")

  if(length(grep("_", lineage_info$lineage)) != 0) {
    if (length(which(lineages$subclade %in% letters)) != 0) {
      lineages<-lineages[-c(which(lineages$subclade %in% letters)),]
    }
  }

  clades<-unique(lineages$subclade)

  if(length(grep("\\.", clades)) != 0 ) {
    clades<-clades[-c(grep("\\.", clades))]
  }


  lineage<-lineage_info$lineage[-c(grep("_", lineage_info$lineage))]
  cols<-RColorBrewer::brewer.pal(9, "Blues")
  pal<-colorRampPalette(c(cols))
  pal<-rev(pal(length(lineage)))
  lineage_info$colour[-c(grep("_", lineage_info$lineage))]<-pal

  for (i in 1:length(clades)) {
    lineage<-grep(clades[i], lineage_info$lineage)
    cols<-RColorBrewer::brewer.pal(3, Colours[i])
    pal<-colorRampPalette(c(cols))
    pal<-rev(pal(length(lineage)))
    lineage_info$colour[(grep(clades[i], lineage_info$lineage))]<-pal
  }
  attach(sequence_data)

  # Plot a nice figure to save

  plot_tree<-ggtree::ggtree(tree, colour = "grey50", ladderize = T) %<+% sequence_data +
    ggtree::geom_tippoint(color="grey50", size=4)+
    ggtree::geom_tippoint(ggplot2::aes(color=lineage), size=3)  +
    ggtree::theme(plot.title = ggplot2::element_text(size = 40, face = "bold"))+
    ggtree::scale_color_manual(values=c(lineage_info$colour)) +
    ggtree::theme(legend.position = "none")
  
  genotype<-data.frame(lineage = sequence_data$lineage)
  rownames(genotype)<-sequence_data$ID
  
  plot_tree<-ggtree::gheatmap(plot_tree, genotype, offset=-0.01, width=.1, font.size=3, color = NA,
                              colnames_angle=-45, hjust=0) +
    ggtree::scale_fill_manual(values=c(lineage_info$colour), name="lineage")+
    ggtree::theme(legend.position = "none")

  return(plot_tree)

}


