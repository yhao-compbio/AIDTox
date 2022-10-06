# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2022
## This script extracts info of compounds and their identified VNN paths from AIDTox interpretation results   


## functions
source("src/functions.R");
library(RColorBrewer);
library(pals);
library(pheatmap);


## 0. Input arguments 
compound_path_file	<- "data/compound_select_gene_comptoxai_tox21_interpret/binding/gamma-epsilon_0.001_0.1/compound_select_gene_comptoxai_binding_tox21-rt-viability-hek293-p1_400_rt_33_ps_5_re_0_st_0_xs_20_al_0.5_ld_1e-04_model.pt_gamma-epsilon_0.001_0.1_path_relevance_pv.tsv";
compound_feature_file	<- "/home/yunhao1/project/comptox_ai/data/compound_select_gene_comptoxai_tox21_data/binding/compound_gene_comptoxai_binding_tox21-rt-viability-hek293-p1_MultiSURF_400_whole_data.tsv";
reactome_file		<- "https://raw.githubusercontent.com/yhao-compbio/ontology/master/downloads/reactome/UniProt2Reactome_All_Levels.txt";
tox21_file		<- "/home/yunhao1/project/tox_data/downloads/tox21/tox21_10k_library_info.tsv";
id_map_file		<- "/home/yunhao1/project/target/downloads/id_map/hgnc_gene_names.tsv";
out_file		<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/compound_select_gene_comptoxai_binding_tox21-rt-viability-hek293-p1_400_rt_33_ps_5_re_0_st_0_xs_20_al_0.5_ld_1e-04_model.pt_gamma-epsilon_0.001_0.1_target_path_relevance_pv.tsv";

## 1. Filter the identified VNN paths of each compound by its ComptoxAI-connected target proteins  
# obtain the target proteins connected to each compound by ComptoxAI 
compound_feature_df <- read.delim(file = compound_feature_file, header = T, row.names = 1);
compound_target <- lapply(rownames(compound_feature_df), function(rcfd){
	rcfd_id <- which(compound_feature_df[rcfd,] == 1);
	rcfd_targets <- colnames(compound_feature_df)[rcfd_id];
	rcfd_targets <- setdiff(rcfd_targets, "assay_outcome");
	return(rcfd_targets);
});
names(compound_target) <- rownames(compound_feature_df);
# obtain the identified VNN path for each compound 
compound_path_df <- read.delim(file = compound_path_file, header = T, sep = "\t");
compound_path1 <- group.vector.by.categories(compound_path_df$cid, compound_path_df$path_id);
compound_path1 <- lapply(compound_path1, unique);
# for each compound, obtain the identified VNN paths that go through the ComptoxAI-connected target proteins 
compound_path <- mapply(function(ncp, cp){
	ncp_targets <- compound_target[[ncp]];
	if(length(ncp_targets) == 0)	return(logical(0))
	else{
		ncp_path_match <- lapply(ncp_targets, function(nt) grep(nt, cp))
		ncp_path_ids <- unique(unlist(ncp_path_match))
		return(cp[ncp_path_ids])
	}
}, names(compound_path1), compound_path1);
# build vectors that indicate one-to-one mapping between compound and its filtered identified VNN paths 
cp_len <- sapply(compound_path, length);
cp_compounds <- unlist(mapply(function(x, y) rep(x, y), names(cp_len), cp_len));
cp_paths <- unlist(compound_path);
# for each filtered identified VNN path, extract three elements along the path 
cp_target_path <- mapply(function(cp){
	cp_s <- strsplit(cp, "_")[[1]];
	# extract target protein of current path
	cs_len <- length(cp_s);
	cp_t <- cp_s[[cs_len]];
	# extract lowest-level pathway of current pathway 
	cp_p1 <- cp_s[[cs_len-1]];
	# extract second lowest-level pathway of current pathway 
	p2_id <- min(c(3, cs_len-1));
	cp_p2 <- cp_s[[p2_id]];
	return(c(cp_t, cp_p1, cp_p2)); 	
}, cp_paths);

## 2. Map target proteins of interest to their gene symbol IDs
# read in uniprot-to-symbol mapping data from HGNC 
id_map <- read.delim(file = id_map_file, header = T, sep = "\t");
# remove gene symbols that cannot be mapped to a UniProt ID 
matched_id <- which(id_map$UniProt.ID.supplied.by.UniProt. != "");
id_map <- id_map[matched_id, ];
# obtain mapped symbol IDs of each Uniprot ID 
id_uni_ids <- lapply(id_map$UniProt.ID.supplied.by.UniProt., function(idu) strsplit(idu, ", ")[[1]]);
id_uni_len <- sapply(id_uni_ids, length);
# build full map of Uniprot ID ~ symbol ID
id_uni_sym <- mapply(function(imas, iul){
        rep(imas, iul)
}, id_map$Approved.symbol, id_uni_len, SIMPLIFY = F);
id_full_map_df <- data.frame(unlist(id_uni_sym), unlist(id_uni_ids));
colnames(id_full_map_df) <- c("Approved.symbol", "UniProt.ID.supplied.by.UniProt.");
# map target proteins of interest to their gene symbol IDs
cp_target_gene <- sapply(cp_target_path[1,], function(ctp1){
	ctp1_id <- which(id_full_map_df[,2] %in% ctp1);
	ctp_gene <- id_full_map_df[ctp1_id,1];
	return(ctp_gene);
});
#

## 3. Map compounds of interest to their names 
# read in Tox21 compound info as data frame, convert compound names to lower case, creat map data frame between compound name PubChem ID 
tox21_df <- read.delim(file = tox21_file, header = T, sep = "\t");
tox21_df$cid <- sapply(tox21_df$PUBCHEM_CID, function(tdpc) paste("CID", tdpc, sep = "_"));
tox21_df$compound_name <- tolower(tox21_df$SAMPLE_NAME);
tox21_df <- tox21_df[, c("cid", "compound_name")];
tox21_df <- unique(tox21_df);
# map compounds of interest to their names 
cp_compound_name <- sapply(cp_compounds, function(cc){
	cc_id <- which(tox21_df$cid %in% cc)[[1]];
	cc_name <- tox21_df$compound_name[[cc_id]];
	return(cc_name);
}); 

## 4. Map pathways of interest to their names 
# read in Reactome annotation as data frame, create mapping between pathway ID and pathway name 
reactome_df <- read.delim(file = reactome_file, header = F, sep = "\t");
reactome_df <- unique(reactome_df[ ,c("V2", "V4")]);
reactome_path_name <- reactome_df$V4;
names(reactome_path_name) <- reactome_df$V2;
# map pathways of interest to their names 
cp_path_name1 <- reactome_path_name[cp_target_path[2,]];
cp_path_name2 <- reactome_path_name[cp_target_path[3,]];

## 5. Output the extract info of compounds and their filtered identified VNN paths  
cp_df <- data.frame(cp_paths, cp_compounds, cp_compound_name, cp_target_path[2,], cp_path_name1, cp_target_path[3,], cp_path_name2, cp_target_path[1,], cp_target_gene);
colnames(cp_df) <- c("path_id", "cid", "compound_name", "reactome_id1", "pathway_name1", "reactome_id2", "pathway_name2", "uniprot_id", "gene_name");
write.table(cp_df, file = out_file, sep = "\t", col.names = T, row.names = F, quote = F);
