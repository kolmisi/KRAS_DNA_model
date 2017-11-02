# merge reduced Faura cell cycle model with KRAS model I build after optimization and manual editing

# files
# kras_model_manual_boolnet.sif
# cellcycle_boolnet_HGNC_noCDC20.sif

# SIF tables: kras_model_manual_boolnet_SIF , cellcycle_boolnet_HGNC_noCDC20_SIF
# common targets
common_targets <- intersect(kras_model_manual_boolnet_SIF$targets, cellcycle_boolnet_HGNC_noCDC20_SIF$targets)
kras_model_manual_boolnet_SIF[kras_model_manual_boolnet_SIF$targets %in% common_targets,]
cellcycle_boolnet_HGNC_noCDC20_SIF[cellcycle_boolnet_HGNC_noCDC20_SIF$targets %in% common_targets,]

# variables of reduced cell cycle model
unique(as.vector(as.matrix(cellcycle_boolnet_HGNC_noCDC20_SIF[,c(1,3)])))
# ~ of KRAS model
unique(as.vector(as.matrix(kras_model_manual_boolnet_SIF[,c(1,3)])))

# "merged_network_HGNC" contains interactions between KRAS model and Traynard model, stored in file 
# merged_network_HGNC <- read.table("merged_network_HGNC.sif", sep = "\t", stringsAsFactors = F, header = T)
sort(unique(as.vector(as.matrix(merged_network_HGNC[,c(1,3)]))))
# rewrite 2-level variables
# merged_network_HGNC[apply(merged_network_HGNC, 1, function(r) any( grepl("_b",r) )),]
merged_network_HGNC$factors <- gsub("_b2","",gsub("_b1","",merged_network_HGNC$factors))
merged_network_HGNC$targets <- gsub("_b2","",gsub("_b1","",merged_network_HGNC$targets))
# remove duplicates
if (sum(duplicated(merged_network_HGNC))>0){
  merged_network_HGNC <- merged_network_HGNC[!duplicated(merged_network_HGNC),]
}

# remove nodes that are not in Faure model
nodes_to_remove <- setdiff( unique(as.vector(as.matrix(merged_network_HGNC[,c(1,3)]))),
                            unique(c(as.vector(as.matrix(cellcycle_boolnet_HGNC_noCDC20_SIF[,c(1,3)])), 
                                     as.vector(as.matrix(kras_model_manual_boolnet_SIF[,c(1,3)])))) )
rownames(merged_network_HGNC) <- c()
if (length(unique(c(which(merged_network_HGNC$factors %in% nodes_to_remove ), which(merged_network_HGNC$factors %in% nodes_to_remove ))))>0) {
  merged_network_HGNC <- merged_network_HGNC[-unique(c(which(merged_network_HGNC$factors %in% nodes_to_remove ), 
                                                       which(merged_network_HGNC$targets %in% nodes_to_remove ))),]
}
merged_network_HGNC <- merged_network_HGNC[order(merged_network_HGNC$targets,merged_network_HGNC$factors),]; rownames(merged_network_HGNC)<-c()

# save
# write.table(merged_network_HGNC, "merged_cellcycle_boolnet_HGNC_noCDC20_kras_model_manual_boolnet.sif", sep = "\t", quote = F, row.names = F)

merged_from_boolnet_models <- rbind(cellcycle_boolnet_HGNC_noCDC20_SIF, kras_model_manual_boolnet_SIF)
merged_from_boolnet_models$interaction <- as.integer(merged_from_boolnet_models$interaction)

# new connections between the 2 networks
inter_model_connections <- anti_join(merged_network_HGNC, merged_from_boolnet_models) # , by = c("factors", "targets")
inter_model_connections <- inter_model_connections[order(inter_model_connections$targets,inter_model_connections$factors),]
rownames(inter_model_connections) <- c()

# connections <- connections_from_pypath( kras_model_manual_reggraph_unpacked[,c(4,2,5)],
#  kras_model_manual_boolnet_SIF_unpacked[,c(4,2,5)], pypath_dataframe)[,1:3]

# kras_model_manual_boolnet_SIF, cellcycle_boolnet_HGNC_noCDC20_SIF

inter_model_connections_unpacked <- connections_from_pypath( cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked,
                                                             kras_model_manual_boolnet_SIF_unpacked, pypath_dataframe)[,1:3]

# manually "unpacked" KRAS model: "kras_model_manual_boolnet_SIF_unpacked.sif"
kras_model_manual_boolnet_SIF_unpacked <- read.table("kras_model_manual_boolnet_SIF_unpacked.sif", sep = "\t", stringsAsFactors=F, header=T)
cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked <- read.table("cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked.sif",sep="\t",stringsAsFactors=F,header=T)

inter_model_connections_unpacked <- inter_model_connections_unpacked[order(inter_model_connections_unpacked$target_hgnc,
                                                                           inter_model_connections_unpacked$source_hgnc),]

write.table(inter_model_connections_unpacked, "inter_model_connections_unpacked.sif", sep = "\t", quote = F, row.names = F)
write.table(inter_model_connections, "inter_model_connections.sif", sep = "\t", quote = F, row.names = F)

# cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked <- cellcycle_boolnet_HGNC_noCDC20_SIF %>% mutate(factors=strsplit(as.character(factors),"_")) %>% 
#   unnest(factors) %>% mutate(targets = strsplit(as.character(targets), "_")) %>% unnest(targets)
# cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked <- cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[,c("factors","interaction","targets")]
# 
# for(counter in 1:nrow(cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked)-1) {
#   if (!grepl("[A-Z]",cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter+1,1]) ) {
#     cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter+1,1] <- paste(gsub("[0-9]$","",cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter,1]),
#       cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter+1,1],sep = "")
#   }
#    if (!grepl("[A-Z]",cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter+1,3]) ) {
#      cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter+1,3] <- paste(gsub("[0-9]$","",cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter,3]),
#          cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[counter+1,3],sep = "")
#  }
# }

# CCND in pypath
pypath_dataframe[grepl("CCND[0-9]",pypath_dataframe$source_hgnc ) | grepl("CCND[0-9]",pypath_dataframe$target_hgnc ), 1:3]
# target
pypath_dataframe[grepl("CCND[0-9]",pypath_dataframe$target_hgnc ), 1:3]


# merge two models, SIFs:
# KRAS: kras_model_manual_boolnet_SIF
# cell-cycle: cellcycle_boolnet_HGNC_noCDC20_SIF
# connections: inter_model_connections
# 
# MODELS
# cellcycle_boolnet_HGNC_noCDC20
# kras_boolnet_model

# opened file "kras_cell_cycle_merged_manual_boolnet.txt" to manually merge 2 models, changed truth condition for "CDK1_CCNB1_2_3"
kras_cell_cycle_merged_manual_boolnet_model <- loadNetwork("kras_cell_cycle_merged_manual_boolnet.txt")

# defining perturbations from MIDAS table
# kras_dna_cnolist <- CNOlist("kras_dna_midas_table_observables.csv"); perturbation_list <- function_create_perturbation_list(kras_dna_cnolist)

# matrix for attractors
kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors <- function_asynchronous_multiple_perturbations(kras_cell_cycle_merged_manual_boolnet_model, perturbation_list, c("CCND","CDH1"), 100)

# order and filter for variables
# FUNCTION
function_create_df_from_attractors_matrix <- function(attractors_matrix,model) {
if (nrow(attractors_matrix)==length(model$genes)) {
  z <- attractors_matrix
  z <- rbind(z[-which(grepl("CCN",rownames(z))),], z[which(grepl("CCN",rownames(z))),])
  z <- rbind(z[-which(!grepl("[A-Z]",rownames(z))),], z[which(!grepl("[A-Z]",rownames(z))),])
  # remove rows that are all 0 or 1
  # z[rowSums(z) != 0 & rowSums(z) != ncol(z),]
  kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df <- melt(z[!(rownames(z) %in% c("KRAS", "DNA_damage")),])
  colnames(kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df) <- c("readout_name","perturbation","value")
  kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df
}
}
# apply function
kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df <- function_create_df_from_attractors_matrix(
  kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors, kras_cell_cycle_merged_manual_boolnet_model)
# matrix of inputs (perturbations)
# perturbation_list
# perturb_list_values <- list(c(0,0,0), c(1,0,0), c(1,-1,0), c(1,0,-1), c(1,-1,-1))
# input_matrix <- as.matrix(as.data.frame(perturb_list_values)); colnames(input_matrix) <- c(); rownames(input_matrix) <- c("KRAS_mutant","CHEK1i","MK2i")

# PLOT steady states
# df_to_plot <- kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df
height_value <- 500; width_value <- 500; breaks = -2:1; col = c("blue","white","red");cexCol_val=1.2; cexRow_val=1.2
function_simple_heatmap(kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df, input_matrix,
                        "kras_cell_cycle_merged_manual_boolnet_model_SS", width_value, height_value, breaks, col,cexCol_val,cexRow_val)


###########
# ADDING new connections
nodes1 <- sort(unique(as.vector(as.matrix(cellcycle_boolnet_HGNC_noCDC20_SIF_unpacked[,c(1,3)]))))
nodes2 <- sort(unique(as.vector(as.matrix(kras_model_manual_boolnet_SIF_unpacked[,c(1,3)])))); nodes2 <- nodes2[!(nodes2 %in% intersect(nodes1,nodes2))]
kras_cell_cycle_merged_manual_boolnet_model_genes_unpacked <- sort(c(nodes1,nodes2))

# z <- pypath_dataframe[intersect( 
#   which(apply(pypath_dataframe[,c(1,3)], 1, function(r) any(r %in% nodes1 ))),
#   which(apply(pypath_dataframe[,c(1,3)], 1, function(r) any(r %in% nodes2 ))) ), c(1:3)]

# MYT1 -> HGNC name: PKMYT1
search_node <- "PKMYT1"
single_node_connections <- pypath_dataframe[intersect( 
  which(apply(pypath_dataframe[,c(1,3)], 1, function(r) any(r %in% kras_cell_cycle_merged_manual_boolnet_model_genes_unpacked ))),
  which(apply(pypath_dataframe[,c(1,3)], 1, function(r) any(r %in% search_node ))) ), c(1:3)]

relevant_old_nodes_unpacked <- as.vector(as.matrix(single_node_connections[,c(1,3)]))[!as.vector(as.matrix(single_node_connections[,c(1,3)])) %in% search_node]

relevant_old_nodes <- kras_cell_cycle_merged_manual_boolnet_model$genes[unlist(sapply(lapply(relevant_old_nodes_unpacked, 
                                  function(r) {grepl(r, kras_cell_cycle_merged_manual_boolnet_model$genes)}), which))]

# add new connections
kras_cell_cycle_merged_manual_boolnet_txt <- read.table("kras_cell_cycle_merged_manual_boolnet.txt",sep=",",stringsAsFactors=F,header=T)
# add sth
kras_cell_cycle_merged_manual_boolnet_txt <- rbind(c("TP53","CDK2_CCNE1_2 & CDK2_CCNA1_2"), kras_cell_cycle_merged_manual_boolnet_txt)
# TP53 -> RB1 connection
kras_cell_cycle_merged_manual_boolnet_txt[kras_cell_cycle_merged_manual_boolnet_txt$targets=="RB1",2] <-
" (!CDK2_CCNA1_2 & !CDK1_CCNB1_2_3 & !CDK4_6_CCND1_2_3 & !CDK2_CCNE1_2) | (CDKN1B & !CDK1_CCNB1_2_3 & !CDK4_6_CCND1_2_3) | (TP53 & !CDK1_CCNB1_2_3)"
filename <- "z.txt"
write.table(kras_cell_cycle_merged_manual_boolnet_txt, filename, sep = ",", quote = F, row.names = F)
kras_cell_cycle_merged_manual_boolnet_model <- loadNetwork(filename)

nodes_init_ON <- c("CCND","CDH1")
# run asynchronous simuls
kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors <- function_asynchronous_multiple_perturbations(
  kras_cell_cycle_merged_manual_boolnet_model, perturbation_list, nodes_init_ON, 100)
# transform to DF
kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df <- function_create_df_from_attractors_matrix(
  kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors, kras_cell_cycle_merged_manual_boolnet_model)

function_simple_heatmap(kras_cell_cycle_merged_manual_boolnet_model_asynchronous_attractors_df, input_matrix,
                "kras_cell_cycle_merged_manual_boolnet_model_SS_TP53", width_value, height_value, breaks, col, cexCol_val, cexRow_val)
