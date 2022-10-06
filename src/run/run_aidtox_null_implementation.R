# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2022
## This script generates shell scripts that run AIDTox model on outcome-shuffled Tox21 datasets under Reactome pathway hierarchy


## functions
source("src/functions.R");


## 0. Input arguments
train_data_folder	<- "/home/yunhao1/project/comptox_ai/data/compound_select_gene_comptoxai_tox21_data_shuffle/binding/";	# folder name of input outcome-shuffled training-testing data files 
test_data_folder	<- "/home/yunhao1/project/comptox_ai/data/compound_select_gene_comptoxai_tox21_data/binding/compound_gene_comptoxai";	# folder name of input validation data files 	
hierarchy_folder	<- "/home/yunhao1/project/ontology/data/reactome/hierarchy_comptoxai_select/";	# folder name of input sorted reactome hierarchy files 
perf_file		<- "data/compound_select_gene_comptoxai_tox21_implementation/compound_select_gene_comptoxai_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv"	# name of optimal model performance file 
outcome_col		<- "assay_outcome";	# name of column that contains Tox21 assay outcome
output_folder		<- "data/compound_select_gene_comptoxai_tox21_null/binding/";	# folder name of output files 
N_cores			<- 20;	# number of CPUs
job_name		<- "aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null";	# job name

## 1. Process input and output files
# obtain the file name and prefix of training/testing datasets
all_train_files <- list.files(train_data_folder);
input_train_id <- sapply(all_train_files, function(atf) length(strsplit(atf, "MultiSURF_400", fixed = T)[[1]]));
all_train_files <- all_train_files[input_train_id == 2];
all_train_datasets <- sapply(all_train_files, function(atf) strsplit(atf, "_", fixed = T)[[1]][[5]]);
unique_train_datasets <- unique(all_train_datasets);
dataset_name1 <- sapply(unique_train_datasets, function(utd) paste("binding", utd, "400", sep = "_")); 
dataset_name2 <- sapply(unique_train_datasets, function(utd) paste("binding", utd, "MultiSURF_400", sep = "_"));
# read in optimal AIDTox performance info of Tox21 datasets as data frame   
perf_df <- read.delim(file = perf_file, header = T, sep = "\t");
pd_id <- which(perf_df$dataset_name %in% dataset_name1);
perf_df <- perf_df[pd_id, ]; 
# obtain hyperparameter setting of optimal AIDTox model 
perf_hp <- mapply(function(pdhs){
	pdhs_s <- strsplit(pdhs, "_")[[1]];
	pdhs_rpr <- paste(pdhs_s[1:6], collapse = "_");
	pdhs_hp <- c(pdhs_rpr, pdhs_s[[4]], pdhs_s[[10]], pdhs_s[[12]], pdhs_s[[14]]);
	return(pdhs_hp);
}, perf_df$hyperparameter_setting);
perf_hp <- t(perf_hp);
# create a folder for each dataset to store returned results of AIDTox model  
data_output_folder <- sapply(dataset_name1, function(dn){
	# name the folder after dataset  
	dn_folder <- paste(output_folder, dn, "/", sep = "");
	system(paste("mkdir", dn_folder, sep = " "));
	return(dn_folder);
});

## 2. Generate commands for jobs 
commands <- mapply(function(dn1, dn2, dof, ph1, ph2, ph3, ph4, ph5){
	# training file
	atf_id <- sapply(all_train_files, function(atf) length(strsplit(atf, dn2)[[1]]));
	dn_train <- all_train_files[atf_id == 2];
	# testing file
	dn_test <- paste(test_data_folder, dn2, "whole_data.tsv_test.tsv", sep = "_");
	# root pathway file 
	ph1_root <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_root.tsv", sep = "");
	# parent/children node connection file 
	ph1_relation <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_knowledge_by_node.tsv", sep = "");
	# node gene number file  
	ph1_size <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_node_size.tsv", sep = "");
	# node layer number file  
	ph1_layer <- paste(hierarchy_folder, dn1, "/", dn1, "_", ph1, "_st_0_layer.tsv", sep = "");
	# put together command that includes files above and learned optimal hyperparameter setting 
	dn_commands <- sapply(dn_train, function(dt){
		dt_train <- paste(train_data_folder, dt, sep = "");
		dt_output <- paste(dof, dt, sep = "");
		dt_command <- paste("python", "src/dtox.py", dt_train, dn_test, outcome_col, ph1_root, ph1_relation, ph1_size, ph1_layer, ph2, dt_output, ph3, ph4, ph5, sep = " ");
		return(dt_command);
	});
	return(dn_commands);
}, dataset_name1, dataset_name2, data_output_folder, perf_hp[,1], perf_hp[,2], perf_hp[,3], perf_hp[,4], perf_hp[,5], SIMPLIFY = F);
commands <- unlist(commands);
# shuffle commands (in order to balance running time of each shell scripts)
ran_id <- sample(1:length(commands), length(commands));
commands <- commands[ran_id];
# write shell scripts for jobs
generate.parallel.bash.files(commands, as.integer(N_cores), job_name, "src/run/");
