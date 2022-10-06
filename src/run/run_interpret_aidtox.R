# !/usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2022
## This script generates shell scripts that runs interpretation procedure on optimal AIDTox models trained for Tox21 datasets


## functions
source("src/functions.R");


## 0. Input arguments
query_data_folder	<- "/home/yunhao1/project/comptox_ai/data/compound_select_gene_comptoxai_tox21_data/binding/compound_gene_comptoxai";	# folder name of input whole Tox21 datasets  
hierarchy_folder	<- "/home/yunhao1/project/ontology/data/reactome/hierarchy_comptoxai_select/";	# folder name of input sorted reactome hierarchy files 
model_folder		<- "data/compound_select_gene_comptoxai_tox21_implementation/binding/"	# folder name of trained DTox models using original Tox21 datasets  
null_folder		<- "data/compound_select_gene_comptoxai_tox21_null/binding/"	# folder name of trained DTox models using outcome-shuffled Tox21 datasets   
perf_file		<- "data/compound_select_gene_comptoxai_tox21_implementation/compound_select_gene_comptoxai_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv";	# name of optimal model performance file 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome
output_folder		<- "data/compound_select_gene_comptoxai_tox21_interpret/binding/";	# folder name of output files 
N_cores			<- 1;	# number of CPUs
job_name		<- "interpret_dtox_compound_select_gene_comptoxai_tox21_implementation";	# job name

## 1. Process input and output files
# list all files in folder that contains optimal AIDTox models using outcome-shuffled Tox21 datasets 
all_null_files <- list.files(null_folder, recursive = TRUE);
# select optimal AIDTox model files (ending with 'model.pt') 
model_id <- sapply(all_null_files, function(anf) length(strsplit(anf, "model")[[1]]));
all_null_files <- all_null_files[model_id == 2];
# obtain the file name and prefix of training/testing datasets
all_model_datasets <- sapply(all_null_files, function(anf) strsplit(anf, "_", fixed = T)[[1]][[2]]);
unique_model_datasets <- unique(all_model_datasets);
dataset_name1 <- sapply(unique_model_datasets, function(umd) paste("binding", umd, "400", sep = "_"));
dataset_name2 <- sapply(unique_model_datasets, function(umd) paste("binding", umd, "MultiSURF_400", sep = "_"));
# match trained AIDTox model files with Tox21 datasets 
null_name_files <- sapply(dataset_name1, function(dn){
	# select files that contain AIDTox model trained for the current dataset    
	anf_id <- sapply(all_null_files, function(anf) length(strsplit(anf, dn)[[1]]));
	dn_files <- all_null_files[anf_id == 2];
	dn_files <- sapply(dn_files, function(dfs) paste(null_folder, dfs, sep = ""));
	# write names of selected files to output file  
	dn_out <- paste(null_folder, dn, "_null_files.txt", sep = "");
	writeLines(dn_files, dn_out);
	return(dn_out);
});

## 2. Obtain optimal hyperparameter setting of optimal AIDTox model for each Tox21 dataset
# read in optimal performance info of Tox21 datasets as data frame   
perf_df <- read.delim(file = perf_file, header = T, sep = "\t");
pd_id <- which(perf_df$dataset_name %in% dataset_name1);
perf_df <- perf_df[pd_id, ];
# obtain hyperparameter setting of trained optimal model 
perf_hp <- mapply(function(pdhs){
	pdhs_s <- strsplit(pdhs, "_")[[1]];
	pdhs_rpr <- paste(pdhs_s[1:6], collapse = "_");
	pdhs_hp <- c(pdhs_rpr, pdhs_s[[4]], pdhs_s[[10]], pdhs_s[[12]], pdhs_s[[14]]);
	return(pdhs_hp);
}, perf_df$hyperparameter_setting);
perf_hp <- t(perf_hp);

## 3. Generate parts of commands 
part <- NULL 
# generate first part of command that includes input dataset/hierarchy/model file info 
part[[1]] <- mapply(function(dn1, dn2, ph1, ph2, ph3, pdhs, nnf){
	# root pathway file  
	ph1_root <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_root.tsv", sep = "");
	# parent/children node connection file
	ph1_relation <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_knowledge_by_node.tsv", sep = "");
	# node gene number file 
	ph1_size <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_node_size.tsv", sep = "");
	# node layer number file  
	ph1_layer <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_layer.tsv", sep = "");
	# node name file 
	ph1_node <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_node.tsv", sep = "");
	# optimal AIDTox model file 
	dn_model <- paste(model_folder, dn1, "/compound_select_gene_comptoxai_", dn1 , "_", pdhs, "_model.pt", sep = "");
	# whole Tox21 dataset file
	dn_data <- paste(query_data_folder, dn2, "whole_data.tsv", sep = "_"); 
	# put together first part command that includes the files above 
	dn_command <- paste("python", "src/interpret_dtox.py", ph1_root, ph1_relation, ph1_size, ph1_layer, ph2, ph3, ph1_node, dn_model, nnf, dn_data, outcome_col, sep = " "); 
	return(dn_command);
}, dataset_name1, dataset_name2, perf_hp[,1], perf_hp[,2], perf_hp[,3], perf_df$hyperparameter_setting, null_name_files);
# generate second part of command that includes LRP rule and hyperparameter values     
part[[2]] <- c("gamma-epsilon 0.001 0.1");
# create a folder for each rule to store retured interpretation results  
part2_folders <- sapply(part[[2]], function(p2){
	# name the folder after LRP rule paired with hyperparameters  
	p2_s <- strsplit(p2, " ")[[1]]
	p2_folder <- paste(output_folder, paste(p2_s, collapse = "_"), "/", sep = "");
	system(paste("mkdir", p2_folder, sep = " "));
	return(p2_folder);
});

## 4. Generate commands for jobs 
commands1 <- generate.all.possible.hyperparameter.combinations(part);
# add output file name to each command 
commands <- sapply(commands1, function(com){
	# obtain AIDTox model file name of current command  
	com_s <- strsplit(com, " ")[[1]];
	N_cs <- length(com_s);
	com_model_s <- strsplit(com_s[[10]], "/", fixed = T)[[1]];
	com_model <- com_model_s[[length(com_model_s)]];
	# obtain LRP rule and its hyperparameters of current command
	com_rule <- paste(com_s[(N_cs-2):N_cs], collapse = "_");
	# name output file after the optimal AIDTox model file paired with LRP rule 
	com_op_file <- paste(output_folder, com_rule, "/", com_model, "_", com_rule, sep = "");
	# add outut file name to end of current command 
	com_new <- paste(com, com_op_file, sep = " ");
	return(com_new);
});
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
