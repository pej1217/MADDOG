tree <- read.tree(file = "iqtree7980/clade_7980.fasta.contree")
metadata <- read.csv(file="iqtree7980/clade_7980.csv")
sequence_data <- read.csv(file="iqtree7980/clade_7980_replace.csv")
lineage_info <- read.csv(file="iqtree7980/clade_7980_lineage_info_replace.csv")
node_data <- read.csv(file="iqtree7980/clade_7980_node_data_replace.csv")
lineage_info_parent <- function(lineage_info, node_data, tree, metadata, sequence_data) {
  tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
  tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)
  previous<-data.frame(assignment = unique(metadata$assignment), parent = "", n_seqs = NA)
  previous$parent[1]<-""
  for (i in 1:length(previous$assignment)) {
    previous$n_seqs[i]<-length(which(metadata$alignment.name == previous$assignment[i]))
  }
  node_data<-node_data[order(node_data$lineage),]

  node_data$parent<-NA
  node_data$parent[1]<-""


  for (i in 2:length(node_data$node)) {
    if (length(which(node_data$node %in% treeio::ancestor(tree, node_data$node[i]))) == 0) {
      node_data$parent[i]<-""
    } else {
      parent<-node_data$lineage[which(node_data$node %in% treeio::ancestor(tree, node_data$node[i]))]
      node_data$parent[i]<-parent[length(parent)]
    }
  }

  lineage_info$parent<-NA

  for (i in 1:length(lineage_info$lineage)) {
    lineage_info$parent[i]<-node_data$parent[which(node_data$lineage == lineage_info$lineage[i])]

  }
  write.csv(lineage_info,file="iqtree7980/clade_7980_lineage_info_parents.csv")
  return(lineage_info)
}
