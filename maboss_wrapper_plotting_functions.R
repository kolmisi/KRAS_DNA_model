# functions for the script in maboss_wrapper_plotting.R

#######################################################################################
#######################################################################################
# FUNCTIONS

# function to create a folder for a model (.bnd file) with different cfg files (KOs/K-ins and initial states)
function_cfg_file_creator <- function(filename_tag_list,initstate_values_list,input_cfg_file, nodes_to_set_initstate, initstate_values, output_nametag,indices_fixed_genes){
  
  file_counter <- 0; cfg_list <- c()
  # generate different filenames INIT STATE variants for diff. mutants (DNA damage)
  for (is_counter in 1:length(filename_tag_list)) {
    output_nametag <- filename_tag_list[is_counter]; initstate_values <- initstate_values_list[[is_counter]]
    # SET INIT STATES
    input_cfg_file_mutants_perturbs <- function_add_initial_states_rewrite(input_cfg_file, nodes_to_set_initstate, initstate_values, output_nametag)
    file_counter <- file_counter+1; cfg_list[file_counter] <- input_cfg_file_mutants_perturbs
    
    for (counter in 1:length(indices_fixed_genes) ) {
      OFF_mutants <- fixed_genes[indices_fixed_genes[[counter]]<0]
      ON_mutants <- fixed_genes[indices_fixed_genes[[counter]]>0]
      if ( length(ON_mutants)==0 ) {filetag <- paste("_fixed",paste0(OFF_mutants, "off",collapse = "_" ), sep="_")} 
      if ( length(OFF_mutants)==0 ) {filetag <- paste("_fixed",paste0(ON_mutants, "on",collapse = "_" ), sep="_")} 
      if ( length(OFF_mutants)>0 && length(ON_mutants)>0 ) {
        filetag <- paste("_fixed", paste0( ON_mutants, "on",collapse = "_"), paste0(OFF_mutants, "off",collapse = "_" ), sep="_")
      }
      # generate mutant CFG
      filename <- function_add_mutant_nodes(input_cfg_file_mutants_perturbs, ON_mutants, OFF_mutants, filetag)
      # store filenames
      file_counter <- file_counter+1; cfg_list[file_counter] <- filename
    }
  }
  # move into directory, create cfg files
  dir_name <- gsub(".cfg","",input_cfg_file); 
  if (dir.exists(dir_name)==F) {
    dir.create(dir_name) 
  } else {
    dir_name <- paste(dir_name, "_1", sep = ""); dir.create(dir_name) 
  }
  file.rename(cfg_list, paste( dir_name, cfg_list, sep = "/") )
  system(paste("cp", input_bnd_file, dir_name))
}

#######################################################################################

function_add_mutant_nodes <- function(input_cfg_fileName, ON_mutants, OFF_mutants, filename_tag) {
  
  file_as_character_string <- readChar(input_cfg_fileName, file.info(input_cfg_fileName)$size)
  output_file_name <- paste( gsub(".cfg","",input_cfg_fileName),  filename_tag, ".cfg", sep="")
  
  if (ON_mutants!="" && length(ON_mutants)>0){
    for (counter in 1:length(ON_mutants)){
      file_as_character_string <- gsub( paste("\\$d_", ON_mutants[counter], "=1;", sep=""),
                                        paste("\\$d_", ON_mutants[counter], "=0;", sep=""), file_as_character_string )
      if (grepl( paste(ON_mutants[counter], ".istate=.*;", sep = ""), file_as_character_string)){
        file_as_character_string <- gsub( paste(ON_mutants[counter] ,"[.]istate=\\d;", sep=""),
                                          paste(ON_mutants[counter] ,".istate=1;", sep=""), file_as_character_string)
      }
    }
  }
  
  if (OFF_mutants!="" && length(OFF_mutants)>0){
    for (counter in 1:length(OFF_mutants)){
      file_as_character_string <- gsub( paste("\\$u_", OFF_mutants[counter], "=1;", sep=""),
                                        paste("\\$u_", OFF_mutants[counter], "=0;", sep=""), file_as_character_string )
      if (grepl( paste(OFF_mutants[counter], ".istate=.*;", sep = ""), file_as_character_string)){
        file_as_character_string <- gsub( paste(OFF_mutants[counter] ,"[.]istate=\\d;", sep=""),
                                          paste(OFF_mutants[counter] ,".istate=0;", sep=""), file_as_character_string)
      }
    }
  }
  fileConn <- file(output_file_name); write(file_as_character_string, fileConn); close(fileConn)
  output_file_name
}

#######################################################################################

function_add_initial_states_rewrite <- function(input_cfg_file, nodes_to_set_initstate, initstate_values, output_nametag){
  # read in cfg file that will be rewritten
  input_cfg_file_char <- unlist(strsplit(readChar(input_cfg_file, file.info(input_cfg_file)$size), "\n"))
  # rows with init states
  init_state_nodes <- input_cfg_file_char[grepl("istate", input_cfg_file_char )]
  # gsub("=.*","", nodes)
  nodes_to_replace <- init_state_nodes[gsub("\\..*","", init_state_nodes) %in% nodes_to_set_initstate]
  # if order of nodes different in .cfg file, need to reorder initstate values accordingly
  if (identical(order(nodes_to_replace), order(nodes_to_set_initstate))==F){
    initstate_values <- initstate_values[match(gsub("\\..*","", nodes_to_replace), nodes_to_set_initstate)]
  }
  
  init_state_nodes[gsub("\\..*", "", init_state_nodes) %in% nodes_to_set_initstate] <- unlist(lapply(1:length(initstate_values), 
                                                    function(r) {gsub("=.*;", paste("=",initstate_values[r],";",sep = ""), nodes_to_replace[r]) } ))
  # replace rows of to be modified initstates with new values
  input_cfg_file_char[grepl("istate", input_cfg_file_char )] <- init_state_nodes
  
  output_file_name_initstates <- paste( gsub(".cfg","",input_cfg_file), output_nametag, ".cfg", sep = "")
  system(paste("cp ", input_cfg_file, output_file_name_initstates))
  # paste("\n\n", paste0(init_nodes, ".istate=", init_values, ";",collapse="\n"), sep = "")
  write(input_cfg_file_char , output_file_name_initstates )
  output_file_name_initstates
}

#######################################################################################

# currently not used in script above
#
# # with this function we can set the init. states of all other states (apart from those we specifically set)
# # to all OFF, ON or random
# function_add_initial_states <- function(input_cfg_file, nodes_to_set_initstate, initstate_values, all_OFF_ON_random, output_nametag){
#   if (any(grepl(all_OFF_ON_random, c("OFF","ON", "on", "off","random")))) {
#     
#     # get list of nodes by finding first and last row with "$" sign
#     input_cfg_file_char <- readChar(input_cfg_file, file.info(input_cfg_file)$size)
#     nodes <- unlist(strsplit(input_cfg_file_char,"\n"))[grepl("\\$",unlist(strsplit(input_cfg_file_char,"\n")) )]
#     nodes <- unique(gsub("\\$u_","",gsub("\\$d_","",gsub("=.*","", nodes))))
#     
#     if (any(grepl(all_OFF_ON_random, c("ON", "on")))  ){
#       init_values <-rep(1, length(nodes))
#       init_values[nodes %in% nodes_to_set_initstate] <- initstate_values
#       init_nodes <- as.character(nodes)
#     }
#     if (any(grepl(all_OFF_ON_random, c("OFF", "off")) ) ){
#       init_values <-rep(0, length(nodes))
#       init_values[nodes %in% nodes_to_set_initstate] <- initstate_values
#       init_nodes <- as.character(nodes)
#     }
#     if (all_OFF_ON_random=="random"){
#       init_values <- initstate_values
#       init_nodes <- nodes_to_set_initstate
#     }
#     
#     output_file_name_initstates <- paste( gsub(".cfg","",input_cfg_file), output_nametag, ".cfg", sep = "")
#     system(paste("cp ", input_cfg_file, output_file_name_initstates))
#     write( paste("\n\n", paste0(init_nodes, ".istate=", init_values, ";",collapse="\n"), sep = ""), output_file_name_initstates, append = T )
#     output_file_name_initstates
#   }
#   else {
#     print("ERROR! specify init states for other nodes: on / off / random")
#   }
# }

#######################################################################################

# saving phenotype plots

plot_dynamics_df <- function(df,foldername,plot_levels,legend_font_size,plot_width,plot_height,line_thickness) {
  
  tidy_df <- melt( df, id=c("Time"))
  
  if (length(plot_levels)<=9){
    require(RColorBrewer)
    myColors <- brewer.pal(length(plot_levels), "Set1")
    names(myColors) <- plot_levels
    colScale <- scale_colour_manual(name = "variable",values = myColors)
  }
  temporal <- ggplot( tidy_df, aes(x=Time,y=value, color=variable)) + 
    geom_line(size=line_thickness) + ylab("Phenotype probablity") +
    theme( legend.position = "right",legend.title = element_blank(), legend.text=element_text(size=legend_font_size) ) + 
    guides(colour = guide_legend(nrow = length(unique(tidy_df$variable)))) + 
    scale_colour_discrete(drop=TRUE, limits = levels(plot_levels)) 
  if (length(plot_levels)<=9){
    temporal <- temporal+ colScale
  }
  ggsave(filename = paste(foldername, paste(deparse(substitute(df)),"_time_evolution.png", sep = ""), sep = "/"), plot=temporal, width = plot_width, height = plot_height)
  # return(temporal)
  
  barplot_data <- as.data.frame(t( df[ nrow(df), which(df[nrow(df), ] > 0)[which(df[nrow(df), ] > 0)>1]]))
  colnames(barplot_data)[1] <- "value"
  barplot_data$ord <- rownames(barplot_data)
  levels(barplot_data$ord ) <- levels(plot_levels)
  # barplot_data
  barplot_elementary_states <- ggplot(barplot_data, aes(x = ord, y = value, label = percent(value), fill=ord ) ) +
    geom_col() +  geom_text(size = font_size/3, position = position_stack(vjust = 0.5)) + coord_flip() + 
    theme(axis.ticks=element_blank(), panel.grid=element_blank(),  # , axis.title= element_blank()
          legend.position = "none",legend.title = element_blank(), text=element_text(size=font_size) ) + ylab("Phenotype probablity") + xlab("") + 
    scale_fill_discrete(drop=TRUE, limits = levels(plot_levels) )
  
  ggsave(filename = paste(foldername, paste(deparse(substitute(df)),"_barplot.png", sep = ""), sep = "/"), 
         plot=barplot_elementary_states, width = plot_width, height = plot_height) 
}

#######################################################################################

multiple_folder_plotter <- function(dir_name, output_folders, plot_width, plot_height, plotting_threshold,line_thickness,font_size ){
  
  setwd(dir_name)
  
  for (counter in 1:length(output_folders) ) {
    output_folder <- output_folders[counter]
    
    output_filename <- paste(output_folder, gsub("\\.","", output_folder), "_probtraj_table.csv", sep = "")
    
    probtraj_table <- read.csv(output_filename, sep = "\t", stringsAsFactors = F,check.names=FALSE)
    if ( all(is.na(probtraj_table[,ncol(probtraj_table) ])) ) {
      probtraj_table <- probtraj_table[,-ncol(probtraj_table) ]
    }
    probtraj_table_only_probs <- probtraj_table[, c(1:4,  seq(5, ncol(probtraj_table), by=2) )]
    colnames(probtraj_table_only_probs) <- gsub("\\]","",gsub("Prob\\[","", colnames(probtraj_table_only_probs) ) )
    
    probtraj_table_only_probs <- probtraj_table_only_probs[,c(1,5:ncol(probtraj_table_only_probs))]

    # display all phenotype trajectories
    probtraj_table_only_probs_TH <- create_plots_probtraj_table(output_folder, plot_width, plot_height, plotting_threshold,line_thickness,font_size)
    
    elementary_states_dynamics <- from_phenotypes_to_state_probabilities( probtraj_table_only_probs )
    
    plot_levels <- factor(sort(colnames(elementary_states_dynamics)[2:ncol(elementary_states_dynamics)]))
    if (any(grepl("nil",as.character(plot_levels)))) {
      plot_levels <- plot_levels[c(plot_levels[-which(plot_levels=="nil")], plot_levels[which(plot_levels=="nil")])]
    }
    # plot_levels <- factor(sort(nodes[nodes_module_entity$entity=="readout"]))
    plot_dynamics_df(elementary_states_dynamics, output_folder, plot_levels, font_size, plot_width, plot_height,line_thickness)
    
  }
  setwd("..")
  # probtraj_table_only_probs_TH
}

#######################################################################################

# create phenotype plots, extract table (df) of probability trajectories

create_plots_probtraj_table <- function(foldername, plot_width, plot_height, plotting_threshold, line_thickness, font_size) {
  
  output_filename <- paste(foldername, gsub("\\.","", foldername), "_probtraj_table.csv", sep = "")
  probtraj_table <- read.csv(output_filename, sep = "\t", stringsAsFactors = F,check.names=FALSE)
  if ( all(is.na(probtraj_table[,ncol(probtraj_table) ])) ) {
    probtraj_table <- probtraj_table[,-ncol(probtraj_table) ]
  }
  probtraj_table_only_probs <- probtraj_table[, c(1:4,  seq(5, ncol(probtraj_table), by=2) )]
  colnames(probtraj_table_only_probs) <- gsub("\\]","",gsub("Prob\\[","", colnames(probtraj_table_only_probs) ) )
  
  probtraj_table_only_probs <- probtraj_table_only_probs[,c(1,5:ncol(probtraj_table_only_probs))]
  
  columns_to_keep <- which( apply( probtraj_table_only_probs > plotting_threshold, 2, sum)>1 &  # has higher values than threshold
                              ( apply( probtraj_table_only_probs, 2, which.max) > 1 | # highest value is not first AND/OR last value not 0
                                  probtraj_table_only_probs[nrow(probtraj_table_only_probs),] > plotting_threshold) ) 
  # probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>plotting_threshold
  
  probtraj_to_plot <- probtraj_table_only_probs[,columns_to_keep]
  
  # probtraj_table_only_probs[,c(1, which(probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>plotting_threshold)[
  # which(probtraj_table_only_probs[dim(probtraj_table_only_probs)[1],]>plotting_threshold)>4])]
  
  # if (ncol(probtraj_to_plot)>10) {
  #   prominent_states <- order(probtraj_to_plot[nrow(probtraj_to_plot),], decreasing = T)[1:10]
  #   prominent_states <- prominent_states[prominent_states!=1]
  #   probtraj_to_plot <- probtraj_to_plot[,c(1,prominent_states)]
  # }

  probtraj_table_only_nonzero_probs_df <- melt( probtraj_to_plot, id=c("Time"))
  
  # plot dynamics
  column_number <- round(sqrt(length(unique(probtraj_table_only_nonzero_probs_df$variable))))
  temporal <- ggplot(probtraj_table_only_nonzero_probs_df, aes(x=Time,y=value, color=variable)) + 
    geom_line(size=line_thickness) + ylab("Phenotype probability") +
    theme( legend.position = "right",legend.title = element_blank(), text=element_text(size = font_size),
           legend.text=element_text(size=font_size-2) ) +
    guides(fill = guide_legend(nrow = length(unique(probtraj_table_only_nonzero_probs_df$variable))))
  ggsave(filename = paste(foldername, "Phenotypes_probability_time_evolution.png", sep = "/"),
         plot=temporal, width = plot_width, height = plot_height)
  ####
  pie_data <- as.data.frame(t( probtraj_to_plot[ nrow(probtraj_to_plot), 2:ncol(probtraj_to_plot)]))
  colnames(pie_data)[1] <- "value"; pie_data$ord <- rownames(pie_data) # gsub("\\..","\\ & ",gsub("Prob.","",rownames(pie_data)))
  barplot_data <- pie_data[order(pie_data$value, decreasing = T),]
  barplot <- ggplot(barplot_data, aes(x = ord, y = value, label = percent(value), fill=ord ) ) +
    geom_bar(stat = "identity") + geom_text(size = font_size/3, position = position_stack(vjust = 0.5)) + coord_flip() + 
    theme( axis.ticks=element_blank(), panel.grid=element_blank(),  # , axis.title= element_blank()
           legend.position = "none",legend.title = element_blank(), text=element_text(size=font_size) ) + ylab("Phenotype probability") + xlab("")
  ggsave(filename=paste(foldername, "Phenotypes_probability_barplot.png", sep = "/"), plot=barplot, width = plot_width, height = plot_height)
  
  # probtraj_table_only_probs[,probtraj_table_only_probs[nrow(probtraj_table_only_probs),]>0]
  probtraj_to_plot
}

#######################################################################################

# this function splits "phenotypes' (simultaneously activated nodes in a sample trajectory) into their elementary states and sums these across simulations
from_phenotypes_to_state_probabilities <- function(df) {
  
  colindices <- lapply( unique(unlist(strsplit(colnames(df)[2:ncol(df)], "--"))), function(x) which(grepl(x, colnames(df))) )
  elementary_final_states <- cbind(df[,"Time"], data.frame(  matrix(0, nrow = nrow(df), ncol = length(colindices)) ))
  colnames(elementary_final_states) <- c("Time",unique(unlist(strsplit(colnames(df)[2:ncol(df)], "--"))))
  colnames(elementary_final_states) <- gsub(">","",gsub("<","",colnames(elementary_final_states)))
  
  for(counter in 1:length(colindices) ) {
    if (length(colindices[[counter]])>1) {
      elementary_final_states[,counter+1] <- rowSums(df[, colindices[[counter]]])
    } else {
      elementary_final_states[,counter+1] <- df[, colindices[[counter]]]
    }
  }
  if (any(grepl("nil",colnames(elementary_final_states)))){
    if ( all(elementary_final_states$nil[round(nrow(elementary_final_states)/2):nrow(elementary_final_states)]<10e-3) ) {
      elementary_final_states <- elementary_final_states[,-which(colnames(elementary_final_states)=="nil")]
    }
  }
  elementary_final_states
}

#######################################################################################
