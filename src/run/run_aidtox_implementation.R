# !/usr/bin/env Rscript 
## created by Yun Hao @MooreLab 2022
## This script generates shell scripts that run AIDTox model (DTox with chemical-gene connections from ComptoxAI) on Tox21 datasets 


## functions
source("src/functions.R");


## 0. Input arguments
Args			<- commandArgs(T);
input_data_folder	<- Args[1];		# folder name of input training-testing/validation data files  
hierarchy_folder	<- Args[2];		# folder name of input sorted reactome hierarchy files  
output_folder		<- Args[3];		# folder name of output files 
N_cores			<- Args[4];		# number of CPUs 
job_name		<- Args[5];		# job name 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome 

## 1. Process Tox21 dataset files 
# list all files in dataset file folder
all_data_files <- list.files(input_data_folder, recursive = T);
# select files that contains whole datasets 
whole_id <- sapply(all_data_files, function(adf) length(strsplit(adf, "test", fixed = T)[[1]]));
all_train_files <- all_data_files[whole_id == 1];
all_train_files <- sapply(all_train_files, function(atf) paste(input_data_folder, atf, sep = ""));
all_test_files <- sapply(all_train_files, function(atf){
	atf_s <- strsplit(atf, "_train.tsv")[[1]];
	atf_test <- paste(atf_s, "_test.tsv", sep = "");
	return(atf_test);
});
# extract dataset info from file name: connection type, dataset name, number of top predictive gene features  
all_file_info <- mapply(function(atf){
	atf_s <- strsplit(atf, "_")[[1]];
	as_type <- atf_s[[10]];
	as_folder <- paste(atf_s[c(10, 11, 16)], collapse = "_");
	return(c(as_type, as_folder));
}, all_train_files);
# create output folder based on dataset info
all_output_folder <- mapply(function(afi1, afi2){
	afi_folder <- paste(output_folder, afi1, "/", afi2, "/", sep = "");
	system(paste("mkdir", afi_folder, sep = " "));
	return(afi_folder);
}, all_file_info[1,], all_file_info[2,]);

## 2. Process sorted Reactome hierarchy files  
all_file_h_folder <- sapply(all_file_info[2,], function(afi2) paste(hierarchy_folder, afi2, "/", sep = ""));
all_h_file_list <- lapply(all_file_h_folder, function(afhf){  
	# list all files in hierarchy file folder 
	afhf_files <- list.files(afhf);
	# obtain unique name prefix for each sorted Reactome hierarchy 
	all_afhf_heads <- sapply(afhf_files, function(af) strsplit(af, "_st")[[1]][[1]]);
	afhf_heads <- unique(all_afhf_heads);
	# add structure feature indicator to name prefix of each sorted Reactome hierarchy  
	afhf_heads <- sapply(afhf_heads, function(ah) paste(ah, "_st_0", sep = ""));
	return(afhf_heads);
});

## 3. Generate parts of commands  
commands <- mapply(function(atrf, atef, aof, afhf, ahfl){
	ahfl_commands <- sapply(ahfl, function(ah){
		# root pathway file 
		ah_root <- paste(afhf, ah, "_root.tsv", sep = "");
		# parent/children node connection file 
		ah_relation <- paste(afhf, ah, "_knowledge_by_node.tsv", sep = "");
		# node gene number file  
		ah_size <- paste(afhf, ah, "_node_size.tsv", sep = "");
		# node layer number file  
		ah_layer <- paste(afhf, ah, "_layer.tsv", sep = "");
		# minimal size of pathways 
		ah_s <- strsplit(ah, "_")[[1]];
		as_id <- which(ah_s %in% "ps") + 1;
		min_path <- ah_s[[as_id]];
		# output file
		ah_output <- paste(aof, "compound_select_gene_comptoxai_", ah, sep = "");
		# put together command
		ah_command <- paste("python", "src/dtox.py", atrf, atef, outcome_col, ah_root, ah_relation, ah_size, ah_layer, min_path, ah_output, 20, 0.5, 0.0001, sep = " ");
		return(ah_command);
	});
	return(ahfl_commands);
}, all_train_files, all_test_files, all_output_folder, all_file_h_folder, all_h_file_list, SIMPLIFY = F);
commands <- unlist(commands);

## 4. Generate commands for jobs 
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
