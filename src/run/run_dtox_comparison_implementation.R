# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022
## This script generates shell scripts that run DTox on Tox21 datasets, which matches the configuration of optimal AIDTox model, for performance comparison


## functions
source("src/functions.R");


## 0. Input arguments
Args			<- commandArgs(T);
input_optimal_file	<- Args[1];		# name of optimal DTox model performance file  
input_data_folder	<- Args[2];		# folder name of input training-testing/validation data files 
hierarchy_folder	<- Args[3];		# folder name of trained DTox model configuration/performance files
output_folder		<- Args[4];		# folder name of output files   
N_cores			<- Args[5];		# number of CPUs
job_name		<- Args[6];		# job name 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome  

## 1. Process input dataset and DTox model configuration files 
# read in optimal AIDTox model performance file (contains Tox21 dataset name, optimal hyperparameter setting for each dataset, and other info)  
optimal_df <- read.delim(file = input_optimal_file, sep = "\t", header = T);
# obtain disk location of each Tox21 dataset included in the AIDTox optimal performance data frame 
input_data_info <- mapply(function(oddn) strsplit(oddn, "_")[[1]], optimal_df$dataset_name);
N_di <- nrow(input_data_info); 
input_data_files <- list.files(input_data_folder, recursive = T);
whole_data_id <- sapply(input_data_files, function(idf) length(strsplit(idf, "whole_data.tsv")[[1]]));
whole_data_files <- input_data_files[whole_data_id == 1];
whole_data_files <- sapply(whole_data_files, function(wdf) paste(input_data_folder, wdf, sep = ""));
# obtain disk location of AIDTox optimal model configuration/performance file for each Tox21 dataset (named after optimal hyperparameter setting in optimal_df)  
match_data_files <- apply(input_data_info, 2, function(idi){
	idi_id <- sapply(whole_data_files, function(wdf){
		id1_match <- sum(sapply(idi, function(id1) length(grep(id1, wdf))));
		return(id1_match);
	});
	wdf_id <- which(idi_id == N_di);
	wdf_file <- whole_data_files[[wdf_id]];
	return(wdf_file);
});
# generate disk location of output file for DTox comparison implementation of each Tox21 dataset  
output_files <- sapply(match_data_files, function(mdf){
	mdf_s <- strsplit(mdf, "/", fixed = T)[[1]];
	ms_len <- length(mdf_s);
	ms_op <- paste(output_folder, mdf_s[[ms_len-1]], "/", mdf_s[[ms_len]], sep = "");
	return(ms_op);
});

## 2. Generate parts of commands for DTox comparison implementation
part <- NULL
# put together first part of DTox command
part[[1]] <- mapply(function(mdf, odhs, of){
	# combined training and testing dataset file
	mdf_train <- paste(mdf, "train.tsv", sep = "_");
	# validation dataset file 
	mdf_test <- paste(mdf, "test.tsv", sep = "_");
	# root pathway file 
	odhs_s <- strsplit(odhs, "_")[[1]];
	re_id <- which(odhs_s %in% "re");
	odhs_s1 <- paste(odhs_s[1:(re_id+1)], collapse = "_");
	hf_root <- paste(hierarchy_folder, odhs_s1, "_st_0_root.tsv", sep = "");
	# parent/children node connection file 
	hf_relation <- paste(hierarchy_folder, odhs_s1, "_st_0_knowledge_by_node.tsv", sep = "");
	# node gene number file  
	hf_size <- paste(hierarchy_folder, odhs_s1, "_st_0_node_size.tsv", sep = "");
	# node layer number file  
	hf_layer <- paste(hierarchy_folder, odhs_s1, "_st_0_layer.tsv", sep = "");
	# minimal size of pathways 
	os_id <- which(odhs_s %in% "ps") + 1;
	min_path <- odhs_s[[os_id]];
	# output file
	mdf_output <- paste(of, odhs_s1, sep = "_");
	# put together first part that includes the files/parameters above
	mdf_command <- paste("python", "src/dtox.py", mdf_train, mdf_test, outcome_col, hf_root, hf_relation, hf_size, hf_layer, min_path, mdf_output, sep = " ");
	return(mdf_command);
}, match_data_files, optimal_df$hyperparameter_setting, output_files);
# put together second part of DTox command that includes maximal size of node modules, coefficient for auxiliary loss, coefficient for L2 regularization 
part[[2]] <- c("20 0.5 0.0001")

## 3. Generate commands for jobs 
commands <- generate.all.possible.hyperparameter.combinations(part);
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
