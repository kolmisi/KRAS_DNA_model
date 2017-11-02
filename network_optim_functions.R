# network optimization functions

rewrite_reggraph_sif <- function(reggraph_filename) {
  inter_table <- read.csv(reggraph_filename, sep = "\t", header = T)
  if (ncol(inter_table)!=3){
    inter_table <- read.csv(reggraph_filename, sep = " ", header = T)
  }
  inter_table[,2] <- as.character(inter_table[,2]); inter_table[inter_table[,2]=="->",2] <- 1; inter_table[inter_table[,2]=="-|",2] <- -1
  write.table( inter_table, paste0(strsplit(reggraph_filename,  "\\.")[[1]][1], ".sif"), quote=F, sep="\t", col.names=F, row.names = F ) 
}

########################################
# CREATE BOOLNET file from result
function_create_boolnet <- function(processed_model, edge_indices, optim_result, boolnet_filename) {
  select_edges <- processed_model$reacID[edge_indices]
  links_table <- cbind(as.data.frame(t(as.data.frame(strsplit(select_edges,"=")))), colMeans(optim_result$stringsTol)[edge_indices] )
  rownames(links_table) <- c(); colnames(links_table) <- c("factors","targets", "frequency")
  links_table$factors <- as.character(links_table$factors); links_table$targets <- as.character(links_table$targets)
  links_table <- links_table[order(links_table$targets),]
  
  boolnet_table <- aggregate(factors ~ targets, data=links_table, FUN=paste, collapse=" | ")
  boolnet_table$factors <- gsub("\\+"," & ", boolnet_table$factors); boolnet_table$factors <- gsub("\\!"," ! ", boolnet_table$factors)
  
  # input nodes to be included
  factors <- unique(unlist(strsplit(gsub(" \\| ", " ",gsub(" ! ", " ",gsub(" & "," ",unique(as.vector(as.matrix(boolnet_table$factors)))))), " ")))
  factors <- factors[(factors=="")==F]
  input_df <- as.data.frame(cbind(setdiff(factors, boolnet_table$targets), setdiff(factors, boolnet_table$targets)))
  colnames(input_df) <- colnames(boolnet_table); boolnet_table <- rbind(input_df, boolnet_table)
  
  write.table( c("targets, factors", paste(boolnet_table$targets, boolnet_table$factors, sep = ", ")), 
               boolnet_filename, quote=F, sep="\n", row.names=F, col.names=F)
  boolnet_table
}

##################################

# create reggraph with header and SOURCE,INTERACTION,TARGET order
create_regular_reggraph <- function(filename) {
  reggraph_table <- read.table(filename, sep = " ", stringsAsFactors = F, header=T)
  if (colnames(reggraph_table)[1]!="source_hgnc" & colnames(reggraph_table)[3]!="target_hgnc" ){
    boolean_cell_cycle_reggraph <- read.table(filename, sep = " ", stringsAsFactors = F, header = F)
    boolean_cell_cycle_reggraph <- boolean_cell_cycle_reggraph[,3:1]
    colnames(boolean_cell_cycle_reggraph) <- c("source_hgnc","interaction_directed_signed","target_hgnc")
    write.table(boolean_cell_cycle_reggraph, filename, sep = " ", row.names=F, quote = F)
  } else {
    print("This is already a regular-format reggraph!")
  }
}

##################################

#####
# create SIF table from links_table (which still has multiple inputs per row)

function_create_SIF <- function(processed_model, edge_indices, optim_result) {
  select_edges <- processed_model$reacID[edge_indices]
  links_table <- cbind(as.data.frame(t(as.data.frame(strsplit(select_edges,"=")))), 
                       colMeans(optim_result$stringsTol)[edge_indices] )
  rownames(links_table) <- c(); colnames(links_table) <- c("factors","targets", "frequency")
  links_table$factors <- as.character(links_table$factors); links_table$targets <- as.character(links_table$targets)
  links_table <- links_table[order(links_table$targets),]
  
  if (any(apply(links_table, 1, function(r) any(grepl("\\+",r)) ))) {
    links_table <- data.frame(factors = unlist(strsplit(links_table$factors, split = "\\+")),
                              targets=rep(links_table$targets, sapply(strsplit(links_table$factors,split="\\+"),length) ),
                              frequency=rep(links_table$frequency, sapply(strsplit(links_table$factors,split="\\+"),length) ) )
  }
  # link signs
  if (any(apply(links_table, 1, function(r) any(grepl("\\!",r)) ))) {
    links_table$type <- ""
    links_table$type[apply(links_table, 1, function(r) any(grepl("\\!",r)) )] <- "inhib"
    links_table$type[apply(links_table, 1, function(r) any(grepl("\\!",r)) )==F] <- "activ"
    links_table$factors <- gsub("\\!","",links_table$factors)
    links_table <- links_table[order(links_table$factors),]
  }
  links_table <- links_table[,c("factors", "type", "targets", "frequency")]
  links_table <- links_table[order(links_table$targets),]
  links_table
}

#######################################

# rewrite BOOLNET table into SIF file
function_boolnet_to_SIF <- function(boolnet_table){
  z <- boolnet_table %>% mutate(factors = strsplit(factors, " & ")) %>% unnest(factors) %>%
    mutate(factors = strsplit(factors, " \\| ")) %>% unnest(factors)
  z$interaction <- 1; z$interaction[grepl("!",z$factors)] <- -1
  z$factors <- gsub("!","",z$factors); z$factors <- gsub("\\)","",z$factors); z$factors <- gsub("\\(","",z$factors);  
  z <- z[!duplicated(z),]
  z[,c(1,3,2)]
}


#########################################

function_rewrite_model_names <- function(df,modelnames_HGNC_mapping) {
  df_goodnames <- df
  modelnames_HGNC_mapping
  
  for (rowcounter in 1:nrow(df)){
    for (r in 1:nrow(modelnames_HGNC_mapping)) {
      if (r==1) {
        y <- gsub(modelnames_HGNC_mapping$modelname[r], modelnames_HGNC_mapping$HGNC[r], df$factors[rowcounter], ignore.case = F )
      } else {
        y <- gsub(modelnames_HGNC_mapping$modelname[r], modelnames_HGNC_mapping$HGNC[r], y , ignore.case = F)
      }
    }
    df_goodnames[rowcounter,"factors"] <- y
  }
  df_goodnames
}


##########################################
# create list of perturbation from CNOlist

function_create_perturbation_list <- function(cnolist_input) {
  # include perturbations (KO, mutation, inhibitions) into BoolNet model
  cues <- getCues(cnolist_input)
  # remove inhibitors if not applied
  cues[rowSums(cues[,colnames(cues) %in% colnames(getInhibitors(kras_dna_cnolist))])==0,
       colnames(cues) %in% colnames(getInhibitors(kras_dna_cnolist))] <- NA; # cues <- apply(cues, 2, as.numeric)
  # in CellNOptR when inhib is on, this is denoted as a 1. But since in BoolNet i directly represent inhibition as setting target gene to 0,
  # this needs to be inverted
  cues[,colnames(cues) %in% colnames(getInhibitors(kras_dna_cnolist))] <- abs(cues[,colnames(cues) %in% colnames(getInhibitors(kras_dna_cnolist))]-1)
  # to convert all elements of a matrix to numberic: apply(cues, 2, as.numeric)
  cue_list <- apply(cues, 1, list)
  lapply(1:length(cue_list), function(r) {cue_list[[r]][[1]][
    lapply(1:length(cue_list), function(r) {sapply(cue_list[[r]], function(x) is.na(x)==F)})[[r]][,1] ]})
  
  
}
############################################
# simple parameterized heatmap for steady states

function_simple_heatmap <- function(data_df,input_matrix,filename_string,width_value,height_value,breaks,col,cexCol_val,cexRow_val) {
  
  # breaks = -1:1; col = c("white","red")
  # if (class(ss_dataframe)!="matrix") { ss_dataframe <- as.matrix(ss_dataframe) }
  # # if (dev.cur()>1) {dev.off()}
  # heatmap.2( ss_dataframe, breaks=breaks, col=col, Colv = F, 
  #       Rowv = F, dendrogram = "none", col,trace="none", sepwidth = c(0.01,0.01), sepcolor="black",
  #       colsep=0:ncol(attr_perturbed_models_df), rowsep=0:nrow(attr_perturbed_models_df), 
  #       lwid=c(0.2,1), lhei=c(0.2,1), margins=c(10,10), cexCol=0.9, cexRow = 0.95, srtCol=45, key = F, main = title_string )
  
  png(filename=paste(filename_string, "png", sep = "."), width=width_value, height=height_value)
  data_matrix <- rbind(input_matrix, acast(data_df, readout_name~perturbation, value.var="value"))
  
  # breaks = -2:1; col = c("blue","white","red")
  heatmap.2( data_matrix, breaks=breaks, col=col, Colv = F, 
             Rowv = F, dendrogram = "none", col,trace="none", sepwidth = c(0.01,0.01), sepcolor="black",
             colsep=0:ncol(data_matrix), rowsep=0:nrow(data_matrix),
             lwid=c(0.2,4), lhei=c(0.2,4), margin=c(5,13), srtCol=45, key = F,cexCol=cexCol_val,cexRow=cexRow_val) # 
  dev.off()
}

####
# unique steady states and their frequency

function_unique_ss_frequency <- function(all_ss_df) {
  reshaped_matrix <- matrix(all_ss_df$value, length(unique(all_ss_df$perturbation))*length(unique(all_ss_df$readout_name)),
                            length(unique(all_ss_df$bstring_index)) )
  # bitstring index
  colnames(reshaped_matrix) <- unique(all_ss_df$bstring_index)
  # rownames: readout_name+perturbation
  rownames(reshaped_matrix) <- paste(as.character(exp_input_tidy_df$readout_name), exp_input_tidy_df$perturbation, sep=",")
  # UNIQUE STEADY STATES strings
  all_unique_ss <- unique.matrix(t(reshaped_matrix))
  # which steady state how many times
  frequencies_values <- unlist(lapply(1:nrow(all_unique_ss), 
                                      function(x) {sum(colSums(reshaped_matrix == all_unique_ss[x,])==nrow(reshaped_matrix))}))
  all_unique_ss <- as.data.frame(all_unique_ss)
  all_unique_ss$frequencies <- round(frequencies_values/length(unique(all_ss_df$bstring_index))*1000)/1000
  all_unique_ss
  
}

####################################
####################################
####################################
# find rows with the same steady states

function_find_bitstrings_with_steadystates <- function( all_ss_df, all_unique_ss ) {
  # all_ss_df <- CNO_results_all_tidy_df
  reshaped_matrix <- t(matrix(all_ss_df$value, length(unique(all_ss_df$perturbation))*length(unique(all_ss_df$readout_name)),
                              length(unique(all_ss_df$bstring_index)) ))
  # all_unique_ss <-CNO_results_all_unique_ss
  
  ll <- lapply(1:nrow(all_unique_ss), function(x) 
  { which(rowSums(sweep(reshaped_matrix,2,as.matrix(all_unique_ss[x,-ncol(all_unique_ss)]),"-"))==0) } )
  
  # careful, some bitstrings are missing so thse indexes need to be applied to the index of bitstrings
  lapply(1:length(ll), function(x) {unique(all_ss_df$bstring_index)[ll[[x]]]} )
  
}

############################################

# rerun ASYNCHONOUS updating n times, take unique steady states

function_boolnet_asynchronous_updating_attractors <- function(model, init_states, n_simul){
  n_simul <- n_simul
  attr_simuls <- matrix(NA, length(model$genes), n_simul); 
  rownames(attr_simuls) <- model$genes
  for (i in 1:n_simul) {
    if (dim(as.data.frame(plotAttractors(getAttractors(model, type="asynchronous", startStates=list(init_states)))))[2] ==1) {
      attr_simuls[,i] <- unlist(plotAttractors(getAttractors(model, type="asynchronous", startStates=list(init_states))))
      # plotAttractors( getAttractors(model, type="asynchronous", startStates=list(init_states)) )$`1`
    } else {
      print("multiple attractors!")
      
    }
  }
  # unique steady states
  # sum(!duplicated(t(attr_simuls)))
  attr_simuls_unique <- attr_simuls[,!duplicated(t(attr_simuls))]
  attr_simuls_unique
  # attr_simuls_unique[order(rownames(attr_simuls)),]
}

# attr_simuls[,!colSums(apply(attr_simuls,2,is.na))==nrow(attr_simuls)]

###########################################
# simulate KRAS model with 5 exp perturbations
# this is SYNCHRONOUS updating
# 

# "function_single_boolnet_simulation" renamed as "function_boolnet_synchronous_simulation"


function_boolnet_synchronous_simulation <- function(sample_boolnet_model){
  # sample_boolnet_model <- loadNetwork(boolnet_filename)
  
  # generate BOOLNET model, perturbations the same
  perturbed_models <- lapply(1:length(perturbation_list), function(r) 
  {fixGenes(sample_boolnet_model, names(perturbation_list[[r]]), perturbation_list[[r]])})
  # calculate attractors
  attr_perturbed_models <- lapply(1:length(perturbed_models), function(r) {
    # set up initial states
    init_states <- as.table(rep(0,length(perturbed_models[[r]]$genes))); names(init_states) <- perturbed_models[[r]]$genes
    init_states[perturbed_models[[r]]$genes %in% names(perturbation_list[[r]])] <- perturbation_list[[r]]
    # CALCULATE attractors
    getAttractors(perturbed_models[[r]], startStates=list(generateState(perturbed_models[[r]], specs=init_states)), 
                  type="asynchronous")
  })
  attr_perturbed_models_list <- lapply( 1:length(attr_perturbed_models), function(r) {plotAttractors(attr_perturbed_models[[r]])})
  
  # attr_perturbed_models_df format the same
  # as tidy dataframe
  if ( sum(lapply(1:length(attr_perturbed_models_list), function(x) {length(attr_perturbed_models_list[[x]]$`1`)})==0)==0 ){
    nice_df <- cbind(rownames(as.data.frame(attr_perturbed_models_list)), as.data.frame(attr_perturbed_models_list))
    attr_BOOLNET_tidy_df <- melt(nice_df, id.vars = colnames(nice_df)[1])
    #  BOOLNET_results_all_tidy_df <- rbind(BOOLNET_results_all_tidy_df, attr_one_bitstring)
  }
  colnames(attr_BOOLNET_tidy_df) <- colnames(exp_input_tidy_df)
  attr_BOOLNET_tidy_df$perturbation <- as.numeric(attr_BOOLNET_tidy_df$perturbation)
  attr_BOOLNET_tidy_df
  
  # attr_BOOLNET_tidy_df_observable <- attr_BOOLNET_tidy_df[attr_BOOLNET_tidy_df$readout_name %in% exp_input_tidy_df$readout_name,]
  # rownames(attr_BOOLNET_tidy_df_observable) <- c()
  # attr_BOOLNET_tidy_df_observable$perturbation <- exp_input_tidy_df$perturbation
  # 
  # attr_BOOLNET_tidy_df_observable
}

###########################################

standard_saver <- function(datafr) {
  write.table(datafr, deparse(substitute(datafr)), sep = "\t",quote = F, row.names = F)
}


############################################
#
# ggplot heatmap
# output_filename <- paste("df_module_activity_ROMA_score_modules_significant_in_all_pathol_", 
#                          sample_groups[sample_counter], "_", cutoff_groups[cutoff_counter], ".pdf", sep = "")
# pdf( paste(results_folder, output_filename, sep = "/") , width=plot_width, height=plot_height) 
# p <- ggplot(df_module_act, aes(dataset, MODULE) ) + geom_tile(aes(fill = mean_expression), colour = "white" ) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=font_size), 
#         axis.text.y = element_text(size=font_size), # element_text(size=font_size) #  element_blank()
#         axis.title = element_text(size = font_size)) + 
#   labs(fill="") + ylab("") + xlab("") + # + ylab("pathways") +
#   geom_rect(data=frames1, size=frame_width, fill=NA, colour="blue", aes(xmin=dataset-0.5, xmax=dataset + 0.5, ymin=MODULE - 0.5, ymax=MODULE + 0.5) ) +
#   geom_rect(data=frames2, size=frame_width, fill=NA, colour="black", aes(xmin=dataset-0.5, xmax=dataset + 0.5, ymin=MODULE - 0.5, ymax=MODULE + 0.5) ) +
#   scale_fill_gradient(low = "green", high = "red")
# if (grepl("reactome", results_folder)) {
#   p <- p + theme(axis.text.y = element_blank() )
# }
# print(p); dev.off()


#######

# extract the cycle from BOOLNET SIMULATION

function_extract_cycle_boolnet <- function(model, init_states){
  
  path_attractor <- getPathToAttractor(network=model, state=init_states)
  n <- nrow(path_attractor)
  init_state <- path_attractor[n,]
  number_runs <- 100; 
  reruns_matrix <- matrix(0,number_runs,ncol(path_attractor)) # data.frame()
  colnames(reruns_matrix) <- model$genes
  for (i in 1:100) { reruns_matrix[i,] <- as.matrix(stateTransition(model, init_state, type="asynchronous") ) }
  
  l <- apply(abs( sweep(path_attractor, 2, as.numeric(path_attractor[nrow(path_attractor),]),"-"))>0, 1, which)
  cycle_first_column <- which( unlist(lapply(1:length(l), function(r) {
    identical( as.numeric(which(abs(colMeans(reruns_matrix) - init_state)>0)), as.numeric(l[[r]]) ) } 
  ) ) )
  
  cyclic_attractor <- t(path_attractor[cycle_first_column:nrow(path_attractor),])
  cyclic_attractor <- cyclic_attractor[order(rownames(cyclic_attractor)),]
  colnames(cyclic_attractor) <- paste("t",1:ncol(cyclic_attractor), sep = "")
  cyclic_attractor
  
}

################################
# CONNECT models thru PYPATH
function(OLD_NETWORK, NEW_MODULE, PYPATH) {
  # in general form
  indices <- union( intersect(
    which(apply(PYPATH[,c(1,3)], 1, function(r) any(r %in% unique(as.vector(as.matrix(OLD_NETWORK[,c(1,3)]))) ) )),
    which( apply(PYPATH[,c(1,3)], 1, function(r) any(r %in% unique(as.vector(as.matrix(NEW_MODULE[,c(1,3)])))[ 
      !( unique(as.vector(as.matrix(NEW_MODULE[,c(1,3)]))) %in% intersect( unique(as.vector(as.matrix(OLD_NETWORK[,c(1,3)]))), 
                                                                           unique(as.vector(as.matrix(NEW_MODULE[,c(1,3)])))) ) ] ) ) ) ),
    intersect(
      which(apply(PYPATH[,c(1,3)], 1, function(r) any(r %in% unique(as.vector(as.matrix(NEW_MODULE[,c(1,3)]))) ) )),
      which( apply(PYPATH[,c(1,3)], 1, function(r) any(r %in% unique(as.vector(as.matrix(OLD_NETWORK[,c(1,3)])))[ 
        !( unique(as.vector(as.matrix(OLD_NETWORK[,c(1,3)]))) %in% intersect( unique(as.vector(as.matrix(OLD_NETWORK[,c(1,3)]))), 
                                                                              unique(as.vector(as.matrix(NEW_MODULE[,c(1,3)])))) ) ]
      ) ) ) )
  )
  
  connections <- PYPATH[ indices, c("source_hgnc","interaction_directed_signed","target_hgnc","database","ref") ]
  connections$source_hgnc_module <- ""; connections$target_hgnc_module <- ""; connections$source_entity <- ""; connections$target_entity <- ""
  connections <- connections[,c("source_hgnc","interaction_directed_signed","target_hgnc","source_hgnc_module",
                                "target_hgnc_module","source_entity","target_entity","database","ref")]
  connections$no_ref <- unlist(lapply(strsplit(connections$ref, ","), length)); connections <- connections[order(connections$source_hgnc_module),]; 
  rownames(connections) <- c()
  connections
}

################################

# run asynchronous updating with several perturbations
init_genes_ON <- c("CCND","CDH1")

function_asynchronous_multiple_perturbations <- function(model, perturbation_list, init_genes_ON, n_simul){
  
  # matrix for attractors
  matrix_asynchronous_attractors <- matrix(NA,length(model$genes), length(perturbation_list))
  rownames(matrix_asynchronous_attractors) <- model$genes
  # apply perturbations
  for (counter in 1:length(perturbation_list)){
    # introduce perturbations
    model <- fixGenes(model, names(perturbation_list[[counter]]), perturbation_list[[counter]])
    # initial conditions
    kras_cell_cycle_merged_init <- rep(0,length(model$genes)); 
    kras_cell_cycle_merged_init <- as.table(rep(0,length(model$genes))); names(kras_cell_cycle_merged_init) <- model$genes
    # consistent with perturbations
    kras_cell_cycle_merged_init[names(kras_cell_cycle_merged_init) %in% names(perturbation_list[[counter]])] <- perturbation_list[[counter]]
    # consistent with cyclic behavior
    kras_cell_cycle_merged_init[unlist(sapply(lapply(init_genes_ON, function(r) { grepl(r, model$genes) }), which))] <- 1
    # simulate hundred times for each perturb., extract unique steady states
    matrix_asynchronous_attractors[,counter] <- function_boolnet_asynchronous_updating_attractors(model,kras_cell_cycle_merged_init,n_simul)
  }
  
  matrix_asynchronous_attractors <- matrix_asynchronous_attractors[
    order(rownames(matrix_asynchronous_attractors)),]
  matrix_asynchronous_attractors
}