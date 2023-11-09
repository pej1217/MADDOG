#' Sequence Lineage Designation
#'
#' This function designates the lineages of a set of sequences given the tree, alignment, ancestral reconstructions
#' and metadata of the sequences. The sequence IDs must match between all of these.
#' The metadata file must have a column called 'ID' which contains the ID's for all the sequences.
#' It also must have a column called 'year' that contains the collection year
#' of each sequence. Additionally, it must have a column called 'country'. This should contain the country of origin of each sequence.
#' If you don't have this information, you can leave the column blank but it needs to exist.
#' If you have more detailed information about the origin of the sequence (e.g. a state) you can put this in the country column instead of the country.
#' It also needs to have an 'assignment' column which contains any known clade assignments.
#' If you don't have this information, you can leave this column blank but it needs to exist.
#'
#' @param tree A phylogenetic tree
#' @param min.support The minimum bootstrap support required
#' @param alignment The alignment in fasta format that corresponds to the tree
#' @param metadata The metadata corresponding to the sequences in the tree, including "ID" "assignment" "country" and "year"
#' @param ancestral The fasta file of the reconstructed ancestral sequences for the tree, generated by Treetime
#' @return Information about the input sequences and their lineage designation
#' @export
seq_designation <- function(tree, min.support, alignment, metadata, ancestral) {
  tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)
  alignment_matrix <- seqinr::as.matrix.alignment(alignment)
  ancestral_matrix <- seqinr::as.matrix.alignment(ancestral)
  sequences <- 10
  max.support<-100
  # Need it as a matrix for later analyses

  `%notin%` <- Negate(`%in%`)

  #############################################
  #            BOOTSTRAP SUPPORT              #
  #############################################
  # Identify nodes with a bootstrap of over 70 (why would the first ~570 nodes be NA?)
  nodes_70 <- which(tree$node.comment > min.support | tree$node.comment == max.support)
  nodes_70 <- nodes_70 + length(tree$tip.label)

  node_data <- data.frame(Node = nodes_70, n_tips = NA)
  # Make a dataframe ready for values to be put in
  # Fill the first column with the numbers of the nodes identified in the previous steps

  for(i in 1:length(nodes_70)) {
    node_data[i,2] <- length(phangorn::Descendants(tree, nodes_70[i], type = "tips")[[1]])
  }
  # For each node identified in the previous step, count the number of tips descended from that node

  nodes_5 <- node_data[(which(node_data[,2]>= sequences)),]
  # Only carry forwards nodes which have more than 5 tips descended from it
  # This has been identified as the definition for a cluster in previous studies

  #############################################
  #            95% COVERAGE WGS               #
  #############################################
  # Make a dataframe ready to fill with info about number of gaps and N bases, and length of the alignment and sequence
  seq_data <- data.frame(ID = alignment$nam, N = NA, "gap" = NA,
                         Length_before = nchar(alignment$seq[[1]]), Length_after = NA)


  for (i in 1:length(alignment$seq)) {
    seq_data$N[i] <- stringr::str_count(alignment$seq[[i]], pattern = 'n')
    seq_data$gap[i] <- stringr::str_count(alignment$seq[[i]], pattern = '-')
    seq_data$Length_after[i] <- (seq_data$Length_before[i] - seq_data$N[i] - seq_data$gap[i])
  }
  # For each sequence, count the number of n bases and gaps
  # Calculate the length after removing these

  nodes_remove <- phangorn::Ancestors(tree,
                                      (which(tree$tip.label
                                             %in% (seq_data$ID[which(seq_data$Length_after < (seq_data$Length_before * 0.95))])
                                      )),
                                      'all')
  # Identify seqs with less than 95% coverage and corresponding to tip numbers
  # List the ancestor nodes for each of these tip numbers

  if (length(nodes_remove)>0) {
    removes <- nodes_remove[[1]]
    for (i in 2:(length(nodes_remove))) {
      removes <- c(removes, nodes_remove[[i]])
    }
    remove_counts <- data.frame(table(removes))
    # Make a table to count the number the removed sequences descended from each node (e.g. for the deeper nodes, all 10 are descended)

    names(remove_counts) <-c('Node', 'freq')
    # Change the names

    remove_counts$Node <- as.integer(levels(remove_counts$Node))
    # Need to change this, or it creates many levels and causes errors

    new_remove <- remove_counts[which(remove_counts[,1] %in% nodes_5[,1]),]; new_remove
    # Not all nodes are included in the nodes_5 data (some are already excluded)
    # Get rid of the nodes not in the nodes_5 data

    nodes_new<-nodes_5

    for (i in new_remove$Node) {
      nodes_new[which(nodes_new == i), 2] <-(nodes_5[which(nodes_5 == i), 2] - (new_remove[which(new_remove == i), 2]))
    }
    # Take away the number of removed tips from the previous total number of tips calculated for each node

    nodes_5 <- nodes_new[(which(nodes_new[,2] >= sequences)),] # Redo this to remove any that now have less than 5, and write over the old nodes_5 so this is updated with the new tip numbers
  }

  #############################################
  #         DIFFERENCE FROM ANCESTOR          #
  #############################################
  seq_data$Year <- NA # Add another column to the seq data ready to fill in dates
  # Add collection year of each sequence to the table (Use latest, as exact collection not always filled in)
  for (i in 1:length(alignment$seq)) {
    seq_data$Year[i] <- metadata$year[which(metadata$ID == seq_data$ID[i])]
  }

  nodes_5$diff <- NA # Add a column in nodes_5 to count the number of nucleotide differences each cluster has from the old seq

  # For each node of interest, find all the tips
  # Make a note of the differences between the oldest seq in the each cluster/lineage and one of the seqs in the lineage
  # Which differences between the old seq and each seq are shared between all the seqs in the lineage
  # E.g. which lineages show one or more shared nucleotides differences from the ancestor
  # Count these differences and add them to the table to be analysed further (may just be n's)

  nodes_reduced <- data.frame(Nodes = (nodes_5$Node - (1+length(tree$tip.label))))

  for (i in 1:length(nodes_5$Node)) {
    cm <- caper::clade.members(nodes_5$Node[i], tree, include.nodes = F, tip.labels = T)
    seq_cm <- which(seq_data$ID %in% cm)
    old <- which(row.names(ancestral_matrix) == paste("NODE_", (sprintf("%07d", nodes_reduced$Nodes[i])), sep=""))

    tips <- which(row.names(ancestral_matrix) %in% cm)
    x <- which(ancestral_matrix[old,] != ancestral_matrix[(tips[1]),])

    for (j in tips[-c(1)]) {
      x <- x[which(x %in% (which(ancestral_matrix[old,] != ancestral_matrix[j,])))]
      print(x)
      nodes_5$diff[i] <- length(x)
    }
  }

  nodes_diff <- nodes_5[(which(nodes_5[,3]!=0)),] # Get rid of the ones with no differences straight away


  #############################################
  #         OVERLAPPING TIPS REMOVAL          #
  #############################################
  # Add a column to nodes_diff and for each node, count how many of the other nodes of interest are descended from it
  nodes_diff$overlaps <- NA
  for (i in 1:length(nodes_diff$Node)) {
    nodes_diff$overlaps[i] <- length(which((phangorn::allDescendants(tree)[[(nodes_diff[i,1])]]) %in% nodes_diff[,1]))
  }

  # Create a data frame for lineage assignments. Add the tip labels, and a column ready to add the lineage they're assigned to
  lineage_assignments <- data.frame(tip = tree$tip.label, cluster = NA)

  # Order the nodes of interest by the number of times they overlap the other nodes of interest (descending)
  nodes_diff <- nodes_diff[order(-nodes_diff$overlaps),]

  # Add a column called cluster and label the clusters
  nodes_diff$cluster <- c(1:(length(nodes_diff$Node)))


  for (i in 1:(length(nodes_diff$Node))) {
    lineage_assignments[which(lineage_assignments[,1] %in% caper::clade.members(nodes_diff[i,1], tree, include.nodes = F, tip.labels = T)), 2] <- nodes_diff[i,5]
  }
  # For each sequence, see if it's a member of a lineage. If yes, put the number of the cluster in it's lineage assignment
  # Do this in order of the node with the most overlaps to the least, to ensure the assignment is at the lowest possible level
  # E.g. if a sequence is in clusters 1-7, it will appear as 7

  summary <- data.frame(cluster = nodes_diff$cluster, count = NA)

  for (i in 1:(length(summary$cluster))) {
    summary$count[i] <- length(which(lineage_assignments$cluster == summary$cluster[i]))
  }
  # Count the number of sequences assigned to each lineage

  if(length(which(summary$count < 2))<0){
    nodes_diff <- nodes_diff[-c(which(nodes_diff$cluster %in% summary$cluster[(which(summary$count < 2))])),]
  }
  # If any lineages have no sequences in them, remove them as an option from the nodes_diff table

  min <- min(summary$count)

  while (min < 2){
    nodes_diff <- nodes_diff[order(-nodes_diff$overlaps),]
    nodes_diff$cluster <-c(1:(length(nodes_diff$Node)))
    lineage_assignments$cluster <- NA
    for (i in c(1:(length(nodes_diff$Node)))) {
      lineage_assignments[which(lineage_assignments[,1] %in% caper::clade.members((nodes_diff[i,1]), tree, include.nodes = F, tip.labels = T)),2]<-nodes_diff[i,5]
    }
    summary <- data.frame(cluster = nodes_diff$cluster, count = NA)

    for (i in 1:(length(summary$cluster))) {
      summary$count[i] <- length(which(lineage_assignments$cluster == summary$cluster[i]))
    }

    min <- min(summary$count)

    if (min >= 2) {
      print("done")
    } else {
      nodes_diff<-nodes_diff[-c(which(nodes_diff$cluster %in% summary$cluster[(which(summary$count < 2))])), ]
    }
  }
  # Repeat the above steps until there are no clusters with 0 sequences left

  issues<-data.frame(node = nodes_diff$Node, n_tips = nodes_diff$n_tips, cluster = nodes_diff$cluster)

  issues<-issues[order(issues$cluster),]

  issues$parent<-NA
  issues$parent[1]<-""

  for (i in 2:length(issues$node)) {
    if (length(which(issues$node %in% treeio::ancestor(tree, issues$node[i]))) == 0) {
      issues$parent[i]<-""
    } else {
      parent<-issues$cluster[which(issues$node %in% treeio::ancestor(tree, issues$node[i]))]
      issues$parent[i]<-parent[length(parent)]
    }
  }

  issues<-issues[rev(order(issues$parent)),]

  issues$number<-NA

  for (i in 1:length(which(issues$parent %in% 1:1000))) {
    issues$number[i]<-
      issues$n_tips[which(issues$cluster == issues$parent[i])] - issues$n_tips[i]
  }

  if(length(which(issues$number <5)) != 0){
    nodes_diff<-nodes_diff[-c(which(nodes_diff$Node %in% issues$node[which(issues$number < 5)])),]
  }

  for (i in 1:(length(nodes_diff$Node))) {
    lineage_assignments[which(lineage_assignments[,1] %in% caper::clade.members(nodes_diff[i,1], tree, include.nodes = F, tip.labels = T)), 2] <- nodes_diff[i,5]
  }

  nodes_diff$numbers<-1:length(nodes_diff$Node)

  for (i in 1:length(nodes_diff$Node)) {
    lineage_assignments$cluster[which(lineage_assignments$cluster == nodes_diff$cluster[i])]<-
      nodes_diff$numbers[i]
  }

  nodes_diff$cluster<-nodes_diff$numbers

  for(i in 1:length(seq_data$ID)){
    seq_data$cluster[i]<-lineage_assignments$cluster[which(lineage_assignments$tip == seq_data$ID[i])]
  }

  sequence_data<-seq_data
  node_data<-nodes_diff
  node_data<-node_data[order(node_data$overlaps, decreasing = T),]
  sequence_data$previous <- NA
  for (i in 1:length(sequence_data$ID)) {
    sequence_data$previous[i]<-
      metadata$assignment[which(metadata$ID == sequence_data$ID[i])]
  }

  previous_assignments<-data.frame(assignment = unique(sequence_data$previous), node = NA)

  node_data$previous<-NA

  for (i in 1:length(node_data$Node)) {
    clades<-unique(sequence_data$previous[
      which(sequence_data$ID %in% tree$tip.label[c(unlist(
        phangorn::Descendants(tree, node_data$Node[i], type = "tips")))])])

    node_data$previous[i]<-
      paste(c(clades), collapse = ", ")

  }

  for (i in 1:length(previous_assignments$assignment)) {
    previous_assignments$node[i]<-which(node_data$previous == previous_assignments$assignment[i])[1]
    previous_assignments$assignment[i]<-previous_assignments$assignment[i]
  }

  possible_names<-data.frame(names = rep(previous_assignments$assignment, 47))
  previous_assignments$assignment<-paste(previous_assignments$assignment, "_A.1", sep = "")

  for (i in 1:length(previous_assignments$assignment)) {
    node_data$cluster[previous_assignments$node[i]]<-previous_assignments$assignment[i]
  }


  if ((length(which(previous_assignments$node == 1))) == 0) {
    node_data$cluster[1]<-"A.1"
  }

  node_data<-node_data[order(node_data$overlaps, decreasing = T), ]
  node_data$test <- NA
  problem_names<-data.frame(letters = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N",
                                        "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z","AA","AB","AC","AD","AE","AF","AG"
                                        ,"AH","AI","AJ","AK","AL","AM","AN","AO","AP","AQ","AR","AS","AT","AU"))
  possible_names<-possible_names[order(possible_names$names),]
  possible_names<-paste(possible_names, problem_names$letters, sep = "_")
  #issuesªí¥Ü¤£¦bnode1ªº«á¥N¸Ì­±
  issues<-which(node_data$Node %notin% ips::descendants(tree, node_data$Node[1], type = "all", ignore.tip = T))
  x<-1
  y<-1
  numbers<-1
  while (length(issues)>y) {
    issues<-issues[-c(1)]
    node_data$cluster[issues[y]]<-paste(node_data$previous[issues[y]], "_", problem_names$letters[x],".1", sep = "")
    numbers<-c(numbers, issues[y])
    nodes<-ips::descendants(tree, node_data$Node[1], type = "all", ignore.tip = T)
    for (i in 2:length(numbers)){
      nodes<-c(nodes, ips::descendants(tree, node_data$Node[numbers[i]], type = "all", ignore.tip = T))
    }
    issues<-which(node_data$Node %notin% nodes)
    y<-y+1
    
    if (length(grep(problem_names$letters[2], node_data$cluster)) == 0) {
      x<-1
    } else {
      x<-x+1
    }
  }

  fix<-grep(",", node_data$cluster)
  while (length(fix) != 0) {
    letter<-problem_names$letters[(length(which(problem_names$letters %in% node_data$cluster))+1)]
    node_data$cluster<-gsub(node_data$cluster[fix], letter, node_data$cluster)
    fix<-grep(",", node_data$cluster)
  }
  
  
  for (i in 1:length(node_data$Node)) {
    test<-which(node_data$Node %in% ips::descendants(tree, node_data$Node[i], type = "all", ignore.tip = T))
    node_data$test[c(test)] <- paste(node_data$cluster[i], ".1", sep = "")
    
    majors<-which(grepl("_", node_data$test))
    node_data$cluster[c(majors)] <- node_data$test[c(majors)]
    
  
    
    duplicates<-unique(node_data$cluster[duplicated(node_data$cluster)])
    problems<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") == 0)]
    duplicates<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") != 0)]

    for (i in 1:length(duplicates)) {
      test<-which(node_data$cluster == duplicates[i])
      test<-test[-c(1)]
      x<-1
      for (j in 1:length(test)) {
        name<-unlist(stringr::str_split(node_data$cluster[test[j]], "\\."))
        name[length(name)]<-x+as.integer(name[length(name)])
        x<-(x+1)
        node_data$cluster[test[j]]<-paste(c(name), collapse='.' )
      }
    }
  }
  
  max=0
  duplicates_title<-NA
  duplicates_title<-duplicates_title[-c(1)]
    for (m in 1:length(node_data$Node)) {
      if(stringr::str_count(node_data$cluster[m],"\\.")>3){#If "." greater than 3, a new letter should be used to replace it.
        problem_title<-substr(stringr::str_split(node_data$cluster[m],"\\_")[[1]][2],1,7) #need to replace ex:"A.1.1.1".1.1
        for(n in 1:length(node_data$Node)){
          if(substr(stringr::str_split(node_data$cluster[n],"\\_")[[1]][2],1,7)==problem_title&&stringr::str_count(node_data$cluster[n],"\\.")>3){
            duplicates_title<-c(duplicates_title,c(n))#which all have A.1.1.1
          }
          
          
          if(is.na(substr(stringr::str_split(node_data$cluster[n],"\\_")[[1]][2],1,1))==0){#If the cluster is not a numerical value, look for the next available letter.
            used_letter<-substr(stringr::str_split(node_data$cluster[n],"\\_")[[1]][2],1,1)
            if(max>25){
              used_letter<-substr(stringr::str_split(node_data$cluster[n],"\\_")[[1]][2],1,2)#AA has already been used
            }
            used_letter<-which(problem_names$letters==used_letter)
            if(used_letter>max&&length(used_letter)!=0){
              max=used_letter
            }
          }
        }
       
        
        for (a in 1:length(duplicates_title)) {#Replace all that need replacement with new letters
          node_data$cluster[duplicates_title[a]]<-gsub(problem_title,problem_names$letters[max+1],node_data$cluster[duplicates_title[a]])
        }
        duplicates_title<-NA
        duplicates_title<-duplicates_title[-c(1)]
      }
    }


  unclassified<-which(!grepl("_", node_data$cluster))
  unclassified<-unclassified[c(-1)]
  for (i in 1:length(node_data$Node)) {
    test<-which(node_data$Node %in% ips::descendants(tree, node_data$Node[i], type = "all", ignore.tip = T))
    node_data$test[c(test)] <- paste(node_data$cluster[i], ".1", sep = "")
    node_data$cluster[unclassified]<-node_data$test[unclassified]
    duplicates<-unique(node_data$cluster[duplicated(node_data$cluster)])
    problems<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") == 0)]
    duplicates<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") != 0)]

    for (i in 1:length(duplicates)) {
      test<-which(node_data$cluster == duplicates[i])
      test<-test[-c(1)]
      x<-1
      for (j in 1:length(test)) {
        name<-unlist(stringr::str_split(node_data$cluster[test[j]], "\\."))
        name[length(name)]<-x+as.integer(name[length(name)])
        x<-(x+1)
        node_data$cluster[test[j]]<-paste(c(name), collapse='.' )
      }
    }
  }

  max=0
  duplicates_title<-NA
  duplicates_title<-duplicates_title[-c(1)]
  for (m in 1:length(node_data$Node)) {
    if(stringr::str_count(node_data$cluster[m],"\\.")>3){
      problem_title <- substr(node_data$cluster[m], 1, 7)
      
      for(n in 1:length(node_data$Node)){
        if(substr(node_data$cluster[n], 1, 7)==problem_title&&stringr::str_count(node_data$cluster[n],"\\.")>3){
          duplicates_title<-c(duplicates_title,c(n))
        }
        node_data$cluster[3] %in% node_data$cluster[unclassified]
        
        if(is.na(substr(node_data$cluster[n],1,1))==0&&node_data$cluster[n] %in% node_data$cluster[unclassified]){
          used_letter<-substr(node_data$cluster[n],1,1)
          if(max>25){
            used_letter<-substr(node_data$cluster[n],1,2)
          }
          used_letter<-which(problem_names$letters==used_letter)
          if(used_letter>max&&length(used_letter)!=0){
            max=used_letter  
          }
        }
      }
      
      
      for (a in 1:length(duplicates_title)) {
        node_data$cluster[duplicates_title[a]]<-gsub(problem_title,problem_names$letters[max+1],node_data$cluster[duplicates_title[a]])
      }
      duplicates_title<-NA
      duplicates_title<-duplicates_title[-c(1)]
    }
  }


  node_data<-node_data[, -c((grep("test", names(node_data))), grep("previous", names(node_data)))]
  for (i in 1:length(node_data$cluster)) {
    sequence_data$cluster[which(sequence_data$cluster == i)] <- node_data$cluster[i]
  }
  sequence_data<-sequence_data[,-c(4)]
  names(sequence_data)<-c("ID", "n_N", "n_gap", "length", "year", "lineage", "previous")

  #if(length(grep("-", sequence_data$lineage)) != 0){
    #clade<-strsplit(sequence_data$lineage[grep("-", sequence_data$lineage)][1], "-")[[1]][1]
    #sequence_data$lineage[-c(grep("-", sequence_data$lineage))]<-
      #paste(clade, sequence_data$lineage[-c(grep("-", sequence_data$lineage))])
  #}else{
    #if(length(grep("_", sequence_data$lineage)) != 0){
    #clade<-strsplit(sequence_data$lineage[grep("_", sequence_data$lineage)][1], " ")[[1]][1]
    #sequence_data$lineage[-c(grep("_", sequence_data$lineage))]<-
      #paste(clade, sequence_data$lineage[-c(grep("_", sequence_data$lineage))])
  #}}

  sequence_data$lineage[grep("NA", sequence_data$lineage)]<-NA
  return(sequence_data)
}
