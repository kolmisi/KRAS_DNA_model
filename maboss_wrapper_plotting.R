# MaBoss_wrapper_plotter_R_scripts
# (author: Mihaly Koltai (mihaly.koltai@curie.fr))

# below are functions to 
# introduce mutations/perturbations in MaBoSS models easily
# define initial states from simulations
# run simulations of multiple models in one go, create separate folders for them
# create plots from multiple MaBoss simulations and save them in their respective folders

# commenting and examples in process...

# prerequirements: 
# install MaBoSS 2.0 from http://maboss.curie.fr
# source Maboss environment: "source MaBoSS.env" on command line
# load all functions in the function file (some of them are interconnected)
source("maboss_wrapper_plotting_functions.R")

# sample code to use functions
# kras_model_conceptual.cfg is a sample cfg file of a model with KRAS and DNA damage components
input_cfg_file <- "kras_model_conceptual.cfg"; input_bnd_file <- gsub(".cfg",".bnd",input_cfg_file)
fixed_genes <- c("CHEK1","MAPKAPK2")
# to what values these genes would be fixed: -1 is KO, +1 knock-in
indices_fixed_genes <- list(c(0,-1),c(-1,0),c(-1,-1))
# initstate lists: the nodes whose initial states are defined
nodes_to_set_initstate <- c("KRAS","DNA_damage")
# the values that in different versions of the model these nodes assume
initstate_values_list <- list(c(0,0), c(0,1), c(1,0));
# define manually the tags to be added to the different model versions, as they are defined in "initstate_values_list". 
filename_tag_list <- c("_wt","_DNAdam","_KRASmut")

# create folder with cfg files and bnd file
function_cfg_file_creator(filename_tag_list,initstate_values_list,input_cfg_file, nodes_to_set_initstate, initstate_values, output_nametag,indices_fixed_genes)

########################################################################
# RUN SIMULATION

# SPECIFY path to MaBoSS.env and MBSS_FormatTable.pl (set it accordingly)
sh_command_vars <- paste("cd ../.. \ncd MaBoSS-env-2.0/ \nsource MaBoSS.env \ncd .. \ncd", 
                         paste(strsplit(getwd(), "/")[[1]][sapply(strsplit(getwd(), "/"), length)], dir_name, sep = "/"))
# generate list of commands
sh_command <- paste(sh_command_vars, paste0(paste("time MBSS_FormatTable.pl", input_bnd_file, cfg_list ), collapse = "\n"), "cd ..", sep = "\n")
# write shell script with list of commands in folder
write(sh_command, paste(dir_name,"run_simul.sh",sep="/"))
# RUN simulations
system(paste("cd", dir_name, "\n", "sh run_simul.sh") )
# system("rm run_simul.sh")

########################################################################
# PLOTTING
# GET the folder that was changed in last x minutes (set the value accordingly)
time_window <- 10; 
# get the folder names with the data we want to plot
output_folders <- system(paste("find ./*/* -type d -mmin -", time_window, " -print", sep = ""), intern=TRUE); output_folders <- gsub("/.*/","/",output_folders)
# PLOTS
# in each folder 4 plots are generated:
# elementary_states_dynamics_barplot.png: dynamics of elementary states (nodes being ON)
# elementary_states_dynamics_time_evolution.png: steady state of elementary states (nodes being ON)
# Phenotypes_probability_barplot.png: dynamics of "phenotypes" (simultaneously activated nodes)
# Phenotypes_probability_time_evolution.png: steady state of "phenotypes" (simultaneously activated nodes)
#
# parameters for plots. Plotting threshold is the minimal relative frequency of a single state (a node being ON) for elementary_states_*.png
font_size <- 7; plot_width <- 6; plot_height <- 4; font_size<-12; plotting_threshold <- 0.001; line_thickness <- 1
multiple_folder_plotter(dir_name, output_folders, plot_width, plot_height, plotting_threshold, line_thickness, font_size)

# plot and access data of an individual simulation (instead of all of one)
# which folder?
index_folder <- 1; plot_folder <- output_folders[index_folder] # paste(gsub(".bnd","",input_bnd_file),gsub("\\.","",output_folders[index_folder]),sep="")
# enter directory of meta-folder
setwd( gsub(".bnd","",input_bnd_file) )
# plot phenotypes and load probability trajectories of selected simulation into a dataframe
probtraj_table_only_probs_TH <- create_plots_probtraj_table(plot_folder, plot_width, plot_height, plotting_threshold,line_thickness,font_size)
setwd("cd ..")

# access elementary states of phenotypes
elementary_states <- from_phenotypes_to_state_probabilities(probtraj_table_only_probs_TH)
