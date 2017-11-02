# best string
best_bitstring <- kras_dna_SIGNORmodel_results$stringsTol[rowSums(t(abs(t(kras_dna_SIGNORmodel_results$stringsTol) - 
                                                                          kras_dna_SIGNORmodel_results$bString)) )==0,]

ss_size <- dim(t(cutAndPlot(kras_dna_cnolist, signor_network_CNO_model_processed,
                        list(kras_dna_SIGNORmodel_results$bString) )$simResults[[1]]$t1))

attr_CNO_results <- t(cutAndPlot(kras_dna_cnolist, signor_network_CNO_model_processed, 
                                 list(kras_dna_SIGNORmodel_results$bString) )$simResults[[1]]$t1)
rownames(attr_CNO_results) <- colnames(getSignals(kras_dna_cnolist)[timepoint+1]$`1`)
attr_CNO_results <- attr_CNO_results[order(rownames(attr_CNO_results)),]

# read in model
# READ in DATA (MIDAS table)
kras_dna_cnolist <- CNOlist(list.files(pattern = "midas"))
# SET UP data-constrained MODEL
signor_network_CNO_model_processed_mod <- preprocessing(kras_dna_cnolist, 
                                                        readSIF("signor_network_reggraph_MK2_CDC25B_link.sif") )
# best_bitstring
kras_dna_SIGNORmodel_results_mod <- gaBinaryT1(kras_dna_cnolist, signor_network_CNO_model_processed_mod, verbose=TRUE)

attr_CNO_results <- t( cutAndPlot(kras_dna_cnolist, signor_network_CNO_model_processed_mod, 
                                 list(kras_dna_SIGNORmodel_results_mod$bString) )$simResults[[1]]$t1 )
rownames(attr_CNO_results) <- colnames(getSignals(kras_dna_cnolist)[timepoint+1]$`1`)
attr_CNO_results <- attr_CNO_results[order(rownames(attr_CNO_results) ),]

# plot steady states

best_bitstring_result_df <- melt(attr_CNO_results)
# colnames(best_bitstring_result_df) <- colnames(exp_input_tidy_df)

# sum of differences
sum(abs(best_bitstring_result_df$value - exp_input_tidy_df$value))


sample_tidy_df <- BOOLNET_results_all_tidy_df_observable_means

df_to_plot <- cbind(exp_input_tidy_df[,1:2], melt(attr_CNO_results[order(rownames(attr_CNO_results) ),]) )
rownames(df_to_plot) <- c(); colnames(df_to_plot) <- colnames(sample_tidy_df)
# concatenate with inputs, recast as matrix, plot (save PNG) as heatmap
height_value <- 500; width_value <- 500
function_simple_heatmap(df_to_plot, "CNOR_mapkapk2_cdc25b_connection", width_value, height_value )

paste(deparse(substitute(BOOLNET_results_all_unique_ss)), ss_index, "freq",
      BOOLNET_results_all_unique_ss$frequencies[ss_index], sep = "_")

kras_dna_cnolist <- CNOlist( list.files(pattern = "midas") )
signor_network_CNO_model_processed <- preprocessing(kras_dna_cnolist, readSIF("signor_network_reggraph.sif") )
# RUN OPTIMIZATION
kras_dna_SIGNORmodel_results <- gaBinaryT1(kras_dna_cnolist, signor_network_CNO_model_processed )
cutAndPlot(kras_dna_cnolist, signor_network_CNO_model_processed, list(kras_dna_SIGNORmodel_results$bString))

# 
###########
boolnet_filename <- "manual_boolnet.txt"
boolnet_table <- function_create_boolnet(signor_network_CNO_model_processed, kras_dna_SIGNORmodel_results$bString==1, 
                                   kras_dna_SIGNORmodel_results, boolnet_filename)
# introduced MK2 -> CDC25B link
boolnet_model <- loadNetwork("manual_boolnet.txt")
attr_BOOLNET_tidy_df <- function_boolnet_synchronous_simulation(boolnet_model)

# take observables
attr_BOOLNET_tidy_df_observable <- attr_BOOLNET_tidy_df[attr_BOOLNET_tidy_df$readout_name %in% exp_input_tidy_df$readout_name,]
rownames(attr_BOOLNET_tidy_df_observable) <- c()
attr_BOOLNET_tidy_df_observable$perturbation <- exp_input_tidy_df$perturbation

# recast as matrix
data_matrix <- rbind(input_matrix, acast(attr_BOOLNET_tidy_df, readout_name~perturbation, value.var="value"))

# df_to_plot <- cbind(exp_input_tidy_df[,1:2], melt(attr_CNO_results[order(rownames(attr_CNO_results) ),]) )
# rownames(df_to_plot) <- c(); colnames(df_to_plot) <- colnames(sample_tidy_df)

# concatenate with inputs, recast as matrix, plot (save PNG) as heatmap
height_value <- 500; width_value <- 500
function_simple_heatmap(attr_BOOLNET_tidy_df, "boolnet_manual_ss", width_value, height_value )
plotModel(signor_network_CNO_model_processed, kras_dna_cnolist, kras_dna_SIGNORmodel_results$bString) 

#####

png(filename="CNOR_network_beststring.png", width=600, height=600 )
plotModel( signor_network_CNO_model_processed, kras_dna_cnolist, bString=kras_dna_SIGNORmodel_results$bString, 
           graphvizParams=list(nodeWidth=5, fontsize=40) );
dev.off()

################################################################
################################################################

plotModel( signor_network_CNO_model_processed, kras_dna_cnolist, bString=kras_dna_SIGNORmodel_results$bString, 
           graphvizParams=list(nodeWidth=5, fontsize=40) ) # , indexIntegr = c(1,2)
# plot with the removed edges not shown

sif_optimized_model <- function_create_SIF(signor_network_CNO_model_processed, truth_condition, kras_dna_SIGNORmodel_results)
sif_optimized_model <- sif_optimized_model[,1:3]; 
sif_optimized_model$type[sif_optimized_model$type %in% "activ"] <- 1
sif_optimized_model$type[sif_optimized_model$type %in% "inhib"] <- -1
write.table( sif_optimized_model, "sif_optimized_model.sif", quote=F, sep="\t", col.names = F, row.names=F)
# manually introduced MAPKAPK2 -> CDC25B
sif_optimized_model <- readSIF( "sif_optimized_model.sif" )

sif_optimized_model_processed <- preprocessing(kras_dna_cnolist, sif_optimized_model  )

# PLOT model
png(filename="CNOR_network_beststring_MK2_CDC25B.png", width=600, height=600 )
plotModel(sif_optimized_model, kras_dna_cnolist, graphvizParams=list(nodeWidth=5, fontsize=40))
dev.off()

################################################################
# corrected cell cycle model (Traynard 2015, booleanized)

sbml_cellcycle_pauline <- loadSBML("cell_cycle/Traynard_Boolean_MamCC_init.sbml")

init_states <- as.table(rep(0,length(sbml_cellcycle_pauline$genes))); names(init_states) <- sbml_cellcycle_pauline$genes
init_states[names(init_states) %in% c("CycD","Rb","Cdh1","p27","Skp2")] <- 1
path_attractor <- getPathToAttractor(network=sbml_cellcycle_pauline, state=init_states)

n <- nrow(path_attractor)
init_state <- path_attractor[n,]
number_runs <- 100; 
reruns_matrix <- matrix(0,number_runs,ncol(path_attractor)) # data.frame()
colnames(reruns_matrix) <- sbml_cellcycle_pauline$genes
for (i in 1:100) { reruns_matrix[i,] <- as.matrix(stateTransition(sbml_cellcycle_pauline, init_state, type="asynchronous") ) }

l <- apply(abs( sweep(path_attractor, 2, as.numeric(path_attractor[nrow(path_attractor),]),"-"))>0, 1, which)
cycle_first_column <- which( unlist(lapply(1:length(l), function(r) {
  identical( as.numeric(which(abs(colMeans(reruns_matrix) - init_state)>0)), as.numeric(l[[r]]) ) } 
) ) )

cyclic_attractor <- t(path_attractor[cycle_first_column:nrow(path_attractor),])
cyclic_attractor <- cyclic_attractor[order(rownames(cyclic_attractor)),]
# rearrange so first column is where:
# cyclins A,B OFF, Cdh1 ON, Cdc20, UbcH10 OFF
t1 <- which(colSums(abs(cyclic_attractor[c("CycE","CycA","CycB","Cdh1","Cdc20"),]-c(0,0,0,1,0)))==0)

cyclic_attractor <- cyclic_attractor[,c(t1:ncol(cyclic_attractor), 1:(t1-1))]
colnames(cyclic_attractor) <- paste("t",1:ncol(cyclic_attractor), sep = "")

# PLOT
plotSequence(sequence=t(cyclic_attractor[rowSums(cyclic_attractor) != 0 & rowSums(cyclic_attractor) != ncol(cyclic_attractor),]) )

# G2 state: CycA&CycB ON, Cdc20 & UbcH10 OFF
plotSequence(sequence=t( cyclic_attractor[c("CycE","CycA","CycB","Cdh1","Cdc20", "UbcH10"),] ))
###

# save as BoolNet file
saveNetwork(sbml_cellcycle_pauline, file="cellcycle_corrected.txt")
# flip reggraph columns as exported from GINSIM
create_regular_reggraph("cell_cycle/boolean_cell_cycle.reggraph")

rewrite_reggraph_sif("cell_cycle/boolean_cell_cycle.reggraph")
# created SIF file

# create SIF file from manually created BOOLNET model
create_regular_reggraph("manual_boolnet.reggraph")
rewrite_reggraph_sif("manual_boolnet.reggraph")

unique(as.vector(boolean_cell_cycle_reggraph[,c(1,3)]))

cellcycle_model_names <- as.data.frame( unique(as.vector(as.matrix(boolean_cell_cycle_reggraph[,c(1,3)]))) )
cellcycle_model_names[,2] <- toupper(cellcycle_model_names[,1]); rownames(cellcycle_model_names) <- c()
colnames(cellcycle_model_names) <- c("model_name","HGNC")
cellcycle_model_names <- cellcycle_model_names[order(cellcycle_model_names$model_name),]

####################
# fill in HGNC names
# E2F has multiple isoforms
cellcycle_model_names[cellcycle_model_names$model_name=="E2F",2] <- paste0(HGNC_uniprot_fullgenome$HGNC_gene_symbol[
  grepl("^E2F",HGNC_uniprot_fullgenome$HGNC_gene_symbol)], collapse = ",")
# p27 -> CDKN1B
cellcycle_model_names$HGNC[grepl("p27", cellcycle_model_names$model_name)] <- "CDKN1B"

# cyclin isoforms
cyc_names <- paste("CCN",substr(as.character(cellcycle_model_names$model_name[grepl("Cyc",cellcycle_model_names$model_name)]),4,4), sep="")
for (r in 1:length(cyc_names)) {
  cellcycle_model_names[cellcycle_model_names$model_name==cellcycle_model_names$model_name[grepl("Cyc",cellcycle_model_names$model_name)][r],2] <- paste0(HGNC_uniprot_fullgenome$HGNC_gene_symbol[
    grepl(paste(cyc_names[r],"[0-9]$",sep = ""),HGNC_uniprot_fullgenome$HGNC_gene_symbol)], collapse = ",")
}
# RB isoforms
# paste0(HGNC_uniprot_fullgenome$HGNC_gene_symbol[grepl( "^RB[0-9]$" ,HGNC_uniprot_fullgenome$HGNC_gene_symbol)], collapse = ",")
cellcycle_model_names$HGNC[10:11] <- "RB1"

# UBE2C
cellcycle_model_names$HGNC[grepl("UbcH10",cellcycle_model_names$model_name)] <- "UBE2C"
####################
# SIFs/reggraphs
# boolean_cell_cycle_reggraph <- read.table("cell_cycle/boolean_cell_cycle.reggraph", sep = "\t", stringsAsFactors = F, header=T)
# kras_model_manual_reggraph <- read.table("manual_boolnet.reggraph", sep = " ", stringsAsFactors = F, header=T)
# # HGNC names only, if multiple isoforms, all are there. for cyclins both 
# kras_model_manual_reggraph$source_hgnc_all <- kras_model_manual_reggraph$source_hgnc
# kras_model_manual_reggraph$target_hgnc_all <- kras_model_manual_reggraph$target_hgnc
# kras_model_manual_reggraph$source_hgnc_all[kras_model_manual_reggraph$source_hgnc=="CDK1_CYCB"] <- paste("CDK1","CCNB1","CCNB2","CCNB3", sep =",")
# kras_model_manual_reggraph$target_hgnc_all[kras_model_manual_reggraph$target_hgnc=="CDK1_CYCB"] <- paste("CDK1","CCNB1","CCNB2","CCNB3", sep =",")
#
# boolean_cell_cycle_reggraph$source_hgnc_all <- boolean_cell_cycle_reggraph$source_hgnc
# boolean_cell_cycle_reggraph$target_hgnc_all <- boolean_cell_cycle_reggraph$target_hgnc
# # map to HGNC with "cellcycle_model_names"
# boolean_cell_cycle_reggraph$source_hgnc_all <- cellcycle_model_names$HGNC[match(boolean_cell_cycle_reggraph$source_hgnc_all,cellcycle_model_names$model_name)]
# #
# boolean_cell_cycle_reggraph$target_hgnc_all <- cellcycle_model_names$HGNC[match(boolean_cell_cycle_reggraph$target_hgnc_all,cellcycle_model_names$model_name)]
#
# standard_saver(boolean_cell_cycle_reggraph); standard_saver(kras_model_manual_reggraph)
# saved in files "kras_model_manual_reggraph" "boolean_cell_cycle_reggraph"
# manually edited "boolean_cell_cycle_reggraph" to include CDKs with cyclins, these will be complex variables
# boolean_cell_cycle_reggraph <- read.table("boolean_cell_cycle_reggraph", sep = "\t", stringsAsFactors=F, header=T)
##############
#

# DIRECT CONNECTIONS BETWEEN 2 subnetworks
# split compound fields (multiple genes separated by comma)
# kras_model_manual_reggraph_unpacked <- kras_model_manual_reggraph %>% mutate(source_hgnc_all = strsplit(as.character(source_hgnc_all), ",")) %>% unnest(source_hgnc_all) %>% 
#   mutate(target_hgnc_all = strsplit(as.character(target_hgnc_all), ",")) %>% unnest(target_hgnc_all)
# boolean_cell_cycle_reggraph_unpacked <- boolean_cell_cycle_reggraph %>% mutate(source_hgnc_all=strsplit(as.character(source_hgnc_all),",")) %>% unnest(source_hgnc_all) %>%
#   mutate(target_hgnc_all = strsplit(as.character(target_hgnc_all), ",")) %>% unnest(target_hgnc_all)
# # connect networks
# connections <- connections_from_pypath( kras_model_manual_reggraph_unpacked[,c(4,2,5)], 
#  boolean_cell_cycle_reggraph_unpacked[,c(4,2,5)], pypath_dataframe)[,1:3]
# # this will have CDK1, CDK2 etc as indiv genes, but we want CDK*_CCN** complexes
# cylin_names <- unique(as.vector(as.matrix(boolean_cell_cycle_reggraph)))[grepl("CDK[0-9]",
#                                           unique(as.vector(as.matrix(boolean_cell_cycle_reggraph))))]
# # duplicate CDK rows and 
# # unique(as.vector(as.matrix(connections)))[grepl("CDK2$", unique(as.vector(as.matrix(connections))))]
# duplicated_connections <- connections[apply(connections, 1, function(r) any( grepl("CDK2",r) )),]
# # replace by one possible CDK2-CCN* complex
# duplicated_connections[duplicated_connections=="CDK2"] <- cylin_names[grepl("CDK2",cylin_names)][1]
# connections[connections=="CDK2"] <- cylin_names[grepl("CDK2",cylin_names)][2]
# # replace CDK1 with complex name
# connections[connections=="CDK1"] <- cylin_names[grepl("CDK1",cylin_names)]
# # merge 
# if (sum(duplicated(rbind(connections, duplicated_connections))) == 0) {
#   connections <- rbind(connections, duplicated_connections)
# }
# 
# # 2 tables: original name and HGNC names
# z_HGNC <- rbind(kras_model_manual_reggraph[,c(4,2,5)], boolean_cell_cycle_reggraph[,c(4,2,5)] ); colnames(z_HGNC) <- colnames(boolean_cell_cycle_reggraph[,1:3])
# merged_network_HGNC <- rbind(z_HGNC, connections )
# z_modelnames <- rbind(kras_model_manual_reggraph[,1:3], boolean_cell_cycle_reggraph[,1:3] ); 
# merged_network_modelnames <- rbind(z_modelnames, connections )
#  
# # save both tables merged joined vertically
# write.table(cbind(merged_network_HGNC, merged_network_modelnames), "merged_kras_cellcycle.reggraph", sep = "\t",quote = F, row.names = F)
# 
# # mapping cyclin_CDK_complexes
# # cylin_names <- unique(as.vector(as.matrix(merged_network_HGNC)))[grepl("CDK[0-9]", unique(as.vector(as.matrix(merged_network_HGNC))))]
# # merged_network_HGNC[apply(merged_network_HGNC, 1, function(r) any(r %in% unique(as.vector(as.matrix(merged_network_HGNC)))[grepl("CDK[0-9]", unique(as.vector(as.matrix(merged_network_HGNC))))] )),]
# 
# # merged_network_HGNC[apply(merged_network_HGNC, 1, function(r) any( grepl(",", unique(as.vector(as.matrix(merged_network_HGNC))))) ),]
# 
# compound_variables <- unique(as.vector(as.matrix(merged_network_HGNC)))[grepl(",", unique(as.vector(as.matrix(merged_network_HGNC))))]
# # compound_variables[sapply(lapply(lapply(strsplit(compound_variables,","), function(r) {gsub('[0-9]$', '', r)} ), unique),length)==1] <- 
# #   lapply(lapply(strsplit(compound_variables,","), function(r) {gsub('[0-9]$','',r)}),unique)[sapply(lapply(lapply(strsplit(compound_variables,","), function(r) {gsub('[0-9]$','',r)} ),unique),length)==1]
# 
# compound_variables_correct_name <- c()
# for (i in 1:length(compound_variables)) {
#     var <- unlist(strsplit(compound_variables[[i]],","))
#     unique_genes <- unique(gsub('[0-9]$','', var ))
#     if (length(unique_genes) >1) {
#     #isoform_numbers <- 
#     indices <- lapply(1:length(lapply(1:length(unique_genes), function(r) {grepl(unique_genes[r], var)})), 
#                       function(x) var[lapply(1:length(unique_genes), function(r) {grepl(unique_genes[r], var)})[[x]]])
#     isoform_numbers <- unlist(lapply(lapply(indices, function(r) {gsub('[A-Z]','',r)}), function(x) paste(x,collapse = "_")))
#     compound_variables_correct_name[i] <- paste(paste(unique_genes, isoform_numbers, sep = ""), collapse="_")
#     } else {
#       isoform_numbers <- gsub(unique_genes,'', var )
#       compound_variables_correct_name[i] <- paste(unique_genes, paste(isoform_numbers,collapse = "_"), sep="_")
#     }
#     # replace comme-separated entries    
#     merged_network_HGNC[merged_network_HGNC==unlist(compound_variables[[i]])] <- compound_variables_correct_name[i] 
# }
# 
# # 2-level variables should be distinguished
# # lines, columns
# multivalued_indices <- cbind(rep(which(apply(merged_network_modelnames, 1, function(r) any( grepl("\\_b",r) ) ) ), 
#     sapply(lapply(which(apply(merged_network_modelnames, 1, function(r) any( grepl("\\_b",r) ) ) ), 
#           function(x) {which(grepl("\\_b", merged_network_modelnames[x,])) } ), length)),
#   unlist(lapply(which(apply(merged_network_modelnames, 1, function(r) any( grepl("\\_b",r) ) ) ), 
#           function(x) {which(grepl("\\_b", merged_network_modelnames[x,])) } ))
# )
# 
# # add value-tag to variable names in merged_network_HGNC
# for (i in 1:nrow(multivalued_indices)) {
#   merged_network_HGNC[multivalued_indices[i,1],multivalued_indices[i,2]] <- paste(merged_network_HGNC[multivalued_indices[i,1],multivalued_indices[i,2]],
#             strsplit(merged_network_modelnames[multivalued_indices[i,1],multivalued_indices[i,2]], "_")[[1]][2],sep = "_")
# }

dim(merged_network_HGNC)
# manual editing, all "E2F*" should be E2F1_3
# partial match
merged_network_HGNC[apply(merged_network_HGNC, 1, function(r) any( grepl("CDKN1B", r) ) ),]
# complete match
merged_network_HGNC[apply(merged_network_HGNC, 1, function(r) any( r %in% "CDKN1B" ) ),]
# merged_network_HGNC[merged_network_HGNC=="CDKN1B"] <- "CDKN1B_b1"
# merged_network_HGNC[merged_network_HGNC=="RB1"] <- "RB1_b1"
# "CDKN1B"

# mapping of modelnames to HGNC names
merged_network_modelnames_HGNC <- cbind(merged_network_modelnames, merged_network_HGNC[,c(1,3)])

modelnames_HGNC_mapping <- cbind( as.vector(as.matrix(merged_network_modelnames_HGNC[,c(1,3)])), 
                                  as.vector(as.matrix(merged_network_modelnames_HGNC[,4:5])) )
modelnames_HGNC_mapping <- modelnames_HGNC_mapping[duplicated(modelnames_HGNC_mapping)==F,]
modelnames_HGNC_mapping <-modelnames_HGNC_mapping[order(modelnames_HGNC_mapping[,1]),]
colnames(modelnames_HGNC_mapping) <- c("modelname","HGNC")
modelnames_HGNC_mapping <- as.data.frame(modelnames_HGNC_mapping); modelnames_HGNC_mapping$modelname <- as.character(modelnames_HGNC_mapping$modelname)
modelnames_HGNC_mapping$HGNC <- as.character(modelnames_HGNC_mapping$HGNC)
# apply(modelnames_HGNC_mapping, 2, function(r) as.character(r))

# read in BOOLNET files, rewrite variable names
# kras_model_manual_boolnet.txt
# boolean_cell_cycle_boolnet.txt
#
# kras_model_manual_boolnet <- loadNetwork("kras_model_manual_boolnet.txt")
boolean_cell_cycle_boolnet_table <- read.table("boolean_cell_cycle_boolnet.txt", sep = ",", stringsAsFactors = F, header = T)
# rewrite TARGETS with canonical names
boolean_cell_cycle_boolnet_table$targets <- modelnames_HGNC_mapping$HGNC[match(boolean_cell_cycle_boolnet_table$targets, 
                                                                               modelnames_HGNC_mapping$modelname)]
# replace modelnames in FACTORS
boolean_cell_cycle_boolnet_table_goodnames <- function_rewrite_model_names(boolean_cell_cycle_boolnet_table, modelnames_HGNC_mapping)

# KRAS BOOLNET TABLE
kras_model_manual_boolnet_table <- read.table("kras_model_manual_boolnet.txt", sep = ",", stringsAsFactors = F, header = T)
# replace TARGETS
kras_model_manual_boolnet_table$targets <- modelnames_HGNC_mapping$HGNC[match(kras_model_manual_boolnet_table$targets, 
                                                                               modelnames_HGNC_mapping$modelname)]
kras_model_manual_boolnet_table_goodnames <- function_rewrite_model_names(kras_model_manual_boolnet_table, modelnames_HGNC_mapping)

write.table(kras_model_manual_boolnet_table_goodnames, "kras_model_manual_boolnet_table.txt", sep = ", ", quote = F, row.names = F)
write.table(boolean_cell_cycle_boolnet_table_goodnames, "boolean_cell_cycle_boolnet_table.txt", sep = ", ", quote = F, row.names = F)

#################################################

boolean_cell_cycle_boolnet_model <- loadNetwork("boolean_cell_cycle_boolnet_table_NONCYCLIC.txt")
kras_model_manual_boolnet_model <- loadNetwork("kras_model_manual_boolnet_table.txt")

# RE-SIMULATE in BOOLNET
init_states <- as.table(rep(0,length(boolean_cell_cycle_boolnet_model$genes))); names(init_states) <- boolean_cell_cycle_boolnet_model$genes
# match(c("CycD","Rb","Cdh1","p27","Skp2"), modelnames_HGNC_mapping$modelname )
init_ON_cycle <- boolean_cell_cycle_boolnet_model$genes[sort(unlist(sapply(lapply(c("CCND","RB1","CDH1","CDKN","SKP2"), 
                                          function(r) { grepl(r, boolean_cell_cycle_boolnet_model$genes) }), which)))]
init_states[names(init_states) %in% init_ON_cycle] <- 1

# CYCLIC BEHAVIOR
cyclic_attractor <- function_extract_cycle_boolnet(boolean_cell_cycle_boolnet_model, init_states)
#
path_to_attractor <- getPathToAttractor(network=boolean_cell_cycle_boolnet_model, state=init_states)
plotSequence(sequence = getPathToAttractor(network=boolean_cell_cycle_boolnet_model, state=init_states))
# unlist(lapply(5:nrow(path_to_attractor), function(n) {sum(abs(path_to_attractor[n,] - path_to_attractor[n-1,]))}))

# how many variables change
# unlist(lapply(2:ncol(cyclic_attractor), function(n) {sum(abs(cyclic_attractor[,n] - cyclic_attractor[,n-1]))}))

# rearrange so first column is where cyclins A,B OFF, Cdh1 ON, Cdc20, UbcH10 OFF
init_genes_t1 <- boolean_cell_cycle_boolnet_model$genes[(unlist(sapply(lapply(c("CCNA","CCNE","CCNB","CDH1","SKP2"), 
                                    function(r) { grepl(r, boolean_cell_cycle_boolnet_model$genes) }), which)))]
t1 <- which(colSums(abs(cyclic_attractor[rownames(cyclic_attractor) %in% init_genes_t1,]-c(0,0,0,1,1)[order(init_genes_t1)]))==0)
cyclic_attractor <- cyclic_attractor[,c(t1:ncol(cyclic_attractor), 1:(t1-1))]
colnames(cyclic_attractor) <- paste("t", 1:ncol(cyclic_attractor), sep = "")

# PLOT
plotSequence(sequence=t(cyclic_attractor[rowSums(cyclic_attractor) != 0 & rowSums(cyclic_attractor) != ncol(cyclic_attractor),]) )

# G2 state: CycA&CycB ON, Cdc20 & UbcH10 OFF
G2_genes <- boolean_cell_cycle_boolnet_model$genes[unlist(sapply(lapply(c("CCNA","CCNE","CCNB","CDH1","UBE"), 
                                    function(r) { grepl(r, boolean_cell_cycle_boolnet_model$genes) }), which)) ]
plotSequence(sequence=t( cyclic_attractor[rownames(cyclic_attractor) %in% G2_genes,] ))
###

c <- cyclic_attractor[rowSums(cyclic_attractor) != 0 & rowSums(cyclic_attractor) != ncol(cyclic_attractor),]
heatmap.2( c, breaks=breaks, col=col, Colv = F, Rowv = F, dendrogram = "none", col,trace="none", sepwidth = c(0.01,0.01), 
           sepcolor="black", colsep=0:ncol(c), rowsep=0:nrow(c), lwid=c(0.2,4), 
           lhei=c(0.2,4), margin=c(5,13), srtCol=45, key = F)

#############################################
# QUIESCENT PHENOTYPE
boolean_cell_cycle_boolnet_model <- loadNetwork("boolean_cell_cycle_boolnet_table.txt")

# initial state for QUIESCENT PHENOTYPE
init_ON_quiescence <- boolean_cell_cycle_boolnet_model$genes[sort(unlist(sapply(lapply(c("CDKN"), 
                                 function(r) { grepl(r, boolean_cell_cycle_boolnet_model$genes) }), which)))]
init_states <- as.table(rep(0,length(boolean_cell_cycle_boolnet_model$genes)));names(init_states) <- boolean_cell_cycle_boolnet_model$genes
init_states[names(init_states) %in% init_ON_quiescence] <- 1

# ATTRACTOR
p <- getPathToAttractor(network=boolean_cell_cycle_boolnet_model, state=init_states)
# p <- p[,colSums(p) != 0 & colSums(p) != nrow(p)]
plotSequence(sequence = p )

c_heatmap <- t(p); colnames(c_heatmap) <- paste("t", 1:ncol(c_heatmap), sep = "")
c_heatmap[rowSums(c_heatmap) != 0 & rowSums(c_heatmap) != ncol(c_heatmap),]
heatmap.2( c_heatmap, breaks=breaks, col=col, Colv = F, Rowv = F, dendrogram = "none", col,trace="none", sepwidth = c(0.01,0.01), 
           sepcolor="black", colsep=0:ncol(c_heatmap), rowsep=0:nrow(c_heatmap), lwid=c(0.2,4), 
           lhei=c(0.2,4), margin=c(5,13), srtCol=45, key = F)

z <- c()
for (i in 1:100){ 
  z[i] <- sum(abs(p[nrow(p),]-stateTransition(boolean_cell_cycle_boolnet_model_CDK1_notfixed, p[nrow(p),], type="asynchronous")))
}

#############################################
# KRAS model
kras_model_manual_boolnet <- loadNetwork("kras_model_manual_boolnet_table.txt")
# set fixed genes
input_KRAS_MK2_CHEK1 <- c(1,0,0)
newNet <- fixGenes(kras_model_manual_boolnet, c("KRAS","MAPKAPK2","CHEK1"), input_KRAS_MK2_CHEK1)
init_states <- as.table(rep(0,length(newNet$genes))); names(init_states) <- newNet$genes
init_states[names(init_states) %in% names(newNet$fixed[newNet$fixed==1])] <- 1
# calculate attractor
getPathToAttractor(network=newNet, state=init_states)
plotSequence(sequence = getPathToAttractor(network=newNet, state=init_states))

attr_BOOLNET_tidy_df <- function_boolnet_synchronous_simulation(kras_model_manual_boolnet)
# recast as matrix
data_matrix <- rbind(input_matrix, acast(attr_BOOLNET_tidy_df, readout_name~perturbation, value.var="value"))

# df_to_plot <- cbind(exp_input_tidy_df[,1:2], melt(attr_CNO_results[order(rownames(attr_CNO_results) ),]) )
# rownames(df_to_plot) <- c(); colnames(df_to_plot) <- colnames(sample_tidy_df)

# concatenate with inputs, recast as matrix, plot (save PNG) as heatmap
height_value <- 500; width_value <- 500
function_simple_heatmap(attr_BOOLNET_tidy_df, "boolnet_manual_ss", width_value, height_value )

#################################################
# merge two models
boolean_cell_cycle_boolnet_table$targets
kras_model_manual_boolnet_table$targets

#####
# this function NOT correct!
# kras_model_manual_boolnet_SIF <- function_boolnet_to_SIF(kras_model_manual_boolnet_table)
# boolean_cell_cycle_boolnet_SIF <- function_boolnet_to_SIF(boolean_cell_cycle_boolnet_table)
####

# reggraphs from Boolnet files
# toSBML(kras_model_manual_boolnet_model, file="kras_model_manual_boolnet_model.sbml")
# toSBML(boolean_cell_cycle_boolnet_model, file="boolean_cell_cycle_boolnet_model.sbml")
#
# SBML -> REGGRAPH in GINSIM
# GINSIM's reggraph
# kras_model_manual_boolnet_SIF <- read.table("kras_model_manual_boolnet_model.reggraph", sep = " ", stringsAsFactors = F)[,3:1]
# kras_model_manual_boolnet_SIF$V2[kras_model_manual_boolnet_SIF$V2=="-|"] <- -1
# kras_model_manual_boolnet_SIF$V2[kras_model_manual_boolnet_SIF$V2=="->"] <- 1
# write.table(kras_model_manual_boolnet_SIF, "kras_model_manual_boolnet.sif", sep = "\t", quote = F, row.names = F)
# #
# boolean_cell_cycle_boolnet_SIF <- read.table("cellcycle_boolnet_HGNC_noCDC20.reggraph", sep = " ", stringsAsFactors = F)[,3:1]
# boolean_cell_cycle_boolnet_SIF$V2[boolean_cell_cycle_boolnet_SIF$V2=="-|"] <- -1
# boolean_cell_cycle_boolnet_SIF$V2[boolean_cell_cycle_boolnet_SIF$V2=="->"] <- 1
# colnames(kras_model_manual_boolnet_SIF) <- c("factors","interaction","targets"); colnames(boolean_cell_cycle_boolnet_SIF) <- c("factors","interaction","targets")
# write.table(boolean_cell_cycle_boolnet_SIF, "boolean_cell_cycle_boolnet.sif", sep = "\t", quote = F, row.names = F)

# SIF tables: kras_model_manual_boolnet_SIF , boolean_cell_cycle_boolnet_SIF
# common targets
common_targets <- intersect(kras_model_manual_boolnet_SIF$targets, boolean_cell_cycle_boolnet_SIF$targets)
kras_model_manual_boolnet_SIF[kras_model_manual_boolnet_SIF$targets %in% common_targets,]
boolean_cell_cycle_boolnet_SIF[boolean_cell_cycle_boolnet_SIF$targets %in% common_targets,]

colnames(merged_network_HGNC) <- colnames(kras_model_manual_boolnet_SIF)
# merged_network_HGNC$interaction_directed_signed[merged_network_HGNC$interaction_directed_signed=="-|"] <- -1
# merged_network_HGNC$interaction_directed_signed[merged_network_HGNC$interaction_directed_signed=="->"] <- 1
# compare to "merged_network_HGNC"

# kras_cell_cycle_boolnet_SIF <- rbind(boolean_cell_cycle_boolnet_SIF, kras_cell_cycle_boolnet_SIF)
kras_cell_cycle_boolnet_SIF <- kras_cell_cycle_boolnet_SIF[order(kras_cell_cycle_boolnet_SIF$targets),]
merged_network_HGNC <- merged_network_HGNC[order(merged_network_HGNC$targets),]

# new connections between the 2 networks
inter_model_connections <- anti_join(merged_network_HGNC,kras_cell_cycle_boolnet_SIF)
inter_model_connections <- inter_model_connections[order(inter_model_connections$targets,inter_model_connections$factors),]
rownames(inter_model_connections) <- c()

######
# Faure cell cycle model 
# I try removing CDC20 and/or UBE2C
# cellcycle_boolnet <- loadNetwork("cellcycle_boolnet_HGNC_noCDC20.txt")
# cellcycle_boolnet_HGNC_NONCYCLIC.txt
# Faure
model <- loadNetwork("cellcycle_boolnet_HGNC_noCDC20.txt")
# Traynard
model <- loadNetwork("boolean_cell_cycle_boolnet_table_NONCYCLIC.txt")

# INIT CONDITIONS
init_states <- as.table(rep(0,length(model$genes))); names(init_states) <- model$genes
init_states_quiescence <- init_states; init_states_cycle <- init_states
# QUIESCENCE: RB1, CDKN1B, CDH1 are ON
init_states_quiescence[unlist(sapply(lapply(c("CDKN","RB","CDH"), function(r) { grepl(r, model$genes) }), which))] <- 1
init_states_cycle[unlist(sapply(lapply(c("CCND","CDH"), function(r) { grepl(r, model$genes) }), which))] <- 1

# ATTRACTOR for CYCD(t=0)=0
# SYNCHRONOUS
p_synchr_attr_quiescent <- as.data.frame(plotAttractors( getAttractors(model,startStates=list(init_states_quiescence)) ))
p_synchr_attr_quiescent <- as.data.frame(p_synchr_attr_quiescent[order(rownames(p_synchr_attr_quiescent)),])
rownames(p_synchr_attr_quiescent)<-model$genes[order(model$genes)]
colnames(p_synchr_attr_quiescent) <- c()
# getPathToAttractor(network=model, state=init_states_quiescence); plotSequence(sequence=p)

# ASYNCHRONOUS (doesnt work if cyclic attr!)
# 1 simulation
# attr <- plotAttractors( getAttractors(model, type="asynchronous", startStates=list(init_states_quiescence)) )$`1`
# attr <- attr[order(rownames(attr)),]
# RERUN n times
p_asynchr_attr_quiescent <- function_boolnet_asynchronous_updating_attractors(model,init_states_quiescence,100);
p_asynchr_attr_quiescent <- as.data.frame(p_asynchr_attr_quiescent[order(names(p_asynchr_attr_quiescent))]); colnames(p_asynchr_attr_quiescent)<-c()
# compare SYNCHR and ASYNCHR.
z <- cbind(p_asynchr_attr_quiescent, p_synchr_attr_quiescent); z[,!duplicated(t(z))]

#######
# ATTRACTOR for CYCD(t=0)=1
# SYNCHRONOUS UPDATING with initial state CYCD=1
# Path
# p_synchr_attr_cycle <- plotSequence(sequence = getPathToAttractor(network=model, state=init_states_cycle))
# if attractor is SIMPLE
p_synchr_attr_cycle <- as.data.frame(plotAttractors(getAttractors(model, startStates=list(init_states_cycle))))
p_synchr_attr_cycle <- as.data.frame(p_synchr_attr_cycle[order(rownames(p_synchr_attr_cycle)),] ); rownames(p_synchr_attr_cycle) <- model$genes[order(model$genes)]
colnames(p_synchr_attr_cycle) <- c()

####
# IF attractor is CYCLIC
cyclic_attractor <- function_extract_cycle_boolnet(model, init_states_cycle)
init_genes_t1 <- model$genes[(unlist(sapply(lapply(c("CCNA","CCNE","CCNB","CCND","CDH1"), function(r) { grepl(r, model$genes) }), which)))]
t1 <- min(which(colSums(abs(cyclic_attractor[rownames(cyclic_attractor) %in% init_genes_t1,]-c(1,0,0,0,1) ) )==0))
if (t1!=1){
  cyclic_attractor <- cyclic_attractor[,c(t1:ncol(cyclic_attractor), 1:(t1-1))]
  colnames(cyclic_attractor)
}
# variables that change
plotSequence(sequence = t(cyclic_attractor[rowSums(cyclic_attractor)!=0 & rowSums(cyclic_attractor)!=ncol(cyclic_attractor),]))
#####

# ASYNCHRONOUS, if attractor is simple
# getAttractors(model, type="asynchronous", startStates=list(init_states_cycle))
p_asynchr_attr_cycle <- function_boolnet_asynchronous_updating_attractors(model, init_states_cycle, 100)
p_asynchr_attr_cycle <- as.data.frame(p_asynchr_attr_cycle[order(names(p_asynchr_attr_cycle))]); colnames(p_asynchr_attr_cycle)<-c()

# SUMMARY
# reduced (Cdc20 removed) cell cycle model of Faure 2006: it produces desired attractors, without cyclic behavior with the initial conditions
# synchronous and asynchronous gives same result
# initial conditions
# for quiescent behavior we need (tested in GINSIM): (CDH1,RB1)=1, (CCND,CCNE,CCNA,CCNB)=0. Other variables can be random
#
# for cycle: (CCND,CDH1)=1, (CCNA,CCNB)=0. Other variables can be random.

############################################################################################################
# KRAS asynchronous rerun

# # SYNCHRONOUS with function "function_boolnet_synchronous_simulation"
# this function contains several glbal variables referring to perturbations
# kras_asynchronous_simuls_5_perturbs <- matrix(NA, length(kras_boolnet_model$genes)*length(perturbation_list), n_simul)
# for (i in 1:n_simul){
#   kras_asynchronous_simuls_5_perturbs[,i] <- function_boolnet_synchronous_simulation(kras_boolnet_model)[,3]
# }
# kras_asynchronous_simuls_5_perturbs[,!duplicated(t(kras_asynchronous_simuls_5_perturbs))]

# read in model
kras_boolnet_model <- loadNetwork("kras_model_manual_boolnet_table.txt")

# simulate KRAS model synchronously, without function
sums <- rep(NA,length(perturbation_list)); 
kras_boolnet_model_synchronous_attractors <- matrix(NA, length(kras_boolnet_model$genes), length(perturbation_list))
rownames(kras_boolnet_model_synchronous_attractors) <- kras_boolnet_model$genes
for (counter in 1:length(perturbation_list)){
  # fix gene values
  # fixGenes(kras_boolnet_model, c("KRAS","CHEK1","MAPKAPK2"), perturbation_list[[counter]])
  # define init conditions
  kras_init <- rep(0,length(kras_boolnet_model$genes)); kras_init <- as.table(rep(0,length(kras_boolnet_model$genes))); names(kras_init) <- kras_boolnet_model$genes
  kras_init[names(kras_init) %in% names(perturbation_list[[counter]])] <- perturbation_list[[counter]]
  
  p <- plotSequence(sequence = getPathToAttractor( fixGenes(kras_boolnet_model,names(perturbation_list[[counter]]), perturbation_list[[counter]]), 
                                                   state=kras_init) )
  sums[counter] <- sum(abs(p[,ncol(p)] - kras_asynchronous_simuls_5_perturbs[((counter-1)*length(kras_boolnet_model$genes)+1):(counter*length(kras_boolnet_model$genes)),1]))
  kras_boolnet_model_synchronous_attractors[,counter] <- p[,ncol(p)]
}

# ASYNCHRONOUS
kras_boolnet_model_asynchronous_attractors <- matrix(NA,length(kras_boolnet_model$genes),length(perturbation_list))
rownames(kras_boolnet_model_asynchronous_attractors) <- kras_boolnet_model$genes
for (counter in 1:length(perturbation_list)){
model <- fixGenes(kras_boolnet_model,names(perturbation_list[[counter]]), perturbation_list[[counter]])
kras_init <- rep(0,length(kras_boolnet_model$genes)); kras_init <- as.table(rep(0,length(kras_boolnet_model$genes))); names(kras_init) <- kras_boolnet_model$genes
kras_init[names(kras_init) %in% names(perturbation_list[[counter]])] <- perturbation_list[[counter]]
# simulate hundred times for each perturb., extract unique steady states
kras_boolnet_model_asynchronous_attractors[,counter] <- function_boolnet_asynchronous_updating_attractors(model,kras_init,100)
}

# difference between ASYNCHR. and SYNCHR.
sum(colSums(abs(kras_boolnet_model_synchronous_attractors - kras_boolnet_model_asynchronous_attractors)))

# SAME!!! (desired result)

#####################################
# connect "cellcycle_boolnet_HGNC_noCDC20" and "kras_boolnet_model"


