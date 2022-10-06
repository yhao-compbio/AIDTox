# !/usr/bin/env Rscript
## created by Yun Hao @MooreLab 2022
## This script uses Sankey diagram to visualize VNN paths identified for HEK293 cytotoxicity outcome, which connect together compounds, gene features, pathway modules and the outcome. 


## functions
library(networkD3);
library(dplyr);
library(htmlwidgets);


## 0. Input arguments
connection_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/compound_select_gene_comptoxai_binding_tox21-rt-viability-hek293-p1_400_rt_33_ps_5_re_0_st_0_xs_20_al_0.5_ld_1e-04_model.pt_gamma-epsilon_0.001_0.1_target_path_relevance_pv.tsv";

## 1. Visualize tubulin-related VNN connections in Sankey diagram 
# specify plotting parameters  
tubu_genes	<- c("TUBA1A", "TUBA4A", "TUBB2A", "TUBB3", "TUBB4A", "TUBB4B");
path_rela_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/tubulin_pathway_relations.tsv";
group_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/tubulin_node_group.tsv";
# extract tubulin-related VNN paths identified for HEK293 cytotoxicity outcome  
connection_df <- read.delim(file = connection_file, header = T, sep = "\t", stringsAsFactors = F);
cd_id <- which(connection_df$gene_name %in% tubu_genes);
tg_connection_df <- connection_df[cd_id, ]; 
# build data frame that contains tubulin-related connections between compounds and gene features 
tubu_gene_ids <- lapply(tubu_genes, function(tg){
	tg_id <- which(tg_connection_df$gene_name %in% tg);
	return(tg_id);        
});
tubu_gene_compounds <- lapply(tubu_gene_ids, function(tgi) unique(tg_connection_df$compound_name[tgi]));
tgc_len <- sapply(tubu_gene_compounds, length);
gene_target <- unlist(mapply(function(x, y) rep(x, y), tubu_genes, tgc_len, SIMPLIFY = F));
compound_source <- unlist(tubu_gene_compounds);
cg_df <- data.frame(compound_source, gene_target);
# build data frame that contains tubulin-related connections between gene features and pathway modules  
tubu_gene_pathways <- lapply(tubu_gene_ids, function(tgi) unique(tg_connection_df$pathway_name1[tgi]));
tgp_len <- sapply(tubu_gene_pathways, length);
gene_source <- unlist(mapply(function(x, y) rep(x, y), tubu_genes, tgp_len, SIMPLIFY = F));
pathway_target <- unlist(tubu_gene_pathways); 
gp_df <- data.frame(gene_source, pathway_target);
colnames(cg_df) <- colnames(gp_df) <- c("source", "target");
tubu_path_df <- rbind(cg_df, gp_df);
# read in data frame that contains tubulin-related connections among pathway modules 
path_rela_df <- read.delim(file = path_rela_file, header = T, sep = "\t", stringsAsFactors = F);
link_df <- rbind(tubu_path_df, path_rela_df); 
link_df$value <- 1;
# build data frame that contains network node and edge info 
all_nodes <- c(unique(compound_source), tubu_genes, unique(pathway_target), unique(path_rela_df$target));
node_df <- data.frame(node = all_nodes);
link_df$source_id <- match(link_df$source, node_df$node) - 1;
link_df$target_id <- match(link_df$target, node_df$node) - 1;
# build data frame that contains node group info (how nodes are colored)
group_df <- read.delim(file = group_file, header = T, sep = "\t", stringsAsFactors = F);
all_node_group <- all_nodes;
names(all_node_group) <- all_nodes;
all_node_group[group_df$node_name] <- group_df$node_group;
node_df <- data.frame(node = all_nodes, group = all_node_group);
# use Sankey diagram to visualize tubulin-related VNN connections
tubu_sk <- sankeyNetwork(
	Links = link_df, 
	Nodes = node_df,
	Source = "source_id",
	Target = "target_id",
	Value = "value",
	NodeID = "node", 
	NodeGroup = "group",
	fontSize = 30, 
	fontFamily = 'Helvetica',
	nodeWidth = 30, 
	nodePadding = 40,
	height = 700,
	width = 2000,
	sinksRight = TRUE
)
tubu_sk

## 2. Visualize CYP-related VNN connections in Sankey diagram  
# specify plotting parameters
cyp_genes	<- c("CYP1A1", "CYP1A2", "CYP2B6", "CYP2C8", "CYP2D6", "CYP3A4", "CYP11A1", "CYP11B1", "CYP19A1", "CYP27A1", "CYP51A1");
path_rela_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/cytochrome_p450_pathway_relations.tsv";
path_short_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/cytochrome_p450_pathway_short_names.tsv";
drug_atc_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/tox21_10k_library_atc_code.tsv";
atc_manual_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/cytochrome_p450_drugs_atc_manual.tsv";
atc_class_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/cytochrome_p450_drugs_atc_class.tsv"
group_file	<- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/cytochrome_p450_node_group.tsv";
# extract CYP-related VNN paths identified for HEK293 cytotoxicity outcome 
cd_id <- which(connection_df$gene_name %in% cyp_genes);
cg_connection_df <- connection_df[cd_id, ]; 
# build data frame that contains CYP-related connections between compounds and gene features 
cyp_gene_ids <- lapply(cyp_genes, function(cg){
	cg_id <- which(cg_connection_df$gene_name %in% cg);
	return(cg_id);
});
cyp_gene_compounds <- lapply(cyp_gene_ids, function(cgi) unique(cg_connection_df$compound_name[cgi]));
cgc_len <- sapply(cyp_gene_compounds, length);
gene_target <- unlist(mapply(function(x, y) rep(x, y), cyp_genes, cgc_len, SIMPLIFY = F));
compound_source <- unlist(cyp_gene_compounds);
cg_df <- data.frame(compound_source, gene_target);
colnames(cg_df) <- c("compound_name", "target_name");
# build data frame that contains mapping bewteen compounds and their ATC class annotations   
drug_atc_df <- read.delim(file = drug_atc_file, header = T, sep = "\t", stringsAsFactors = F);
cd_atc_df <- merge(drug_atc_df, cg_df, by = "compound_name");
cd_atc3_df <- unique(cd_atc_df[, c("compound_name", "atc_id3")]);
atc_manual_df <- read.delim(file = atc_manual_file, header = T, sep = "\t", stringsAsFactors = F);
cad_id <- which(cd_atc3_df$compound_name %in% atc_manual_df$compound_name);
oad_id <- setdiff(1:nrow(cd_atc3_df), cad_id);
cd_atc_df1 <- rbind(cd_atc3_df[oad_id, ], atc_manual_df);
# build data frame that contains CYP-related connections between ATC classes and gene features  
atc_class_df <- read.delim(file = atc_class_file, header = T, sep = "\t", stringsAsFactors = F);
cd_atc_class_df <- merge(cd_atc_df1, atc_class_df, by = "atc_id3");
cg_class_df <- merge(cg_df, cd_atc_class_df, by = "compound_name");
cg_class_df1 <- unique(cg_class_df[, c("atc_id3_class", "target_name")]);
colnames(cg_class_df1) <- c("source", "target");
cg_class_df1$value <- mapply(function(ccdaic, ccdtn){
	ccdaic_id <- which(cg_class_df$atc_id3_class %in% ccdaic);
	ccdtn_id <- which(cg_class_df$target_name %in% ccdtn);
	ccd_inter_len <- length(intersect(ccdaic_id, ccdtn_id));
	return(ccd_inter_len);
}, cg_class_df1$source, cg_class_df1$target);
# build data frame that contains CYP-related connections between gene features and pathway modules  
cyp_gene_pathways <- lapply(cyp_gene_ids, function(cgi) unique(cg_connection_df$pathway_name2[cgi]));
cgp_len <- sapply(cyp_gene_pathways, length);
gene_source <- unlist(mapply(function(x, y) rep(x, y), cyp_genes, cgp_len, SIMPLIFY = F));
pathway_target <- unlist(cyp_gene_pathways);
gp_df <- data.frame(gene_source, pathway_target);
colnames(gp_df) <- c("source", "target");
gp_df$value <- 1
cyp_path_df <- rbind(cg_class_df1, gp_df);
# read in data frame that contains CYP-related connections among pathway modules  
path_rela_df <- read.delim(file = path_rela_file, header = T, sep = "\t", stringsAsFactors = F);
path_rela_df$value <- 1
link_df <- rbind(cyp_path_df, path_rela_df);
# build data frame that contains network node and edge info  
all_nodes <- c(unique(cg_class_df1$source), cyp_genes, unique(pathway_target), unique(path_rela_df$target));
link_df$source_id <- match(link_df$source, all_nodes) - 1;
link_df$target_id <- match(link_df$target, all_nodes) - 1;
# build data frame that contains node name annotations 
node_short_df <- read.delim(file = path_short_file, header = T, sep = "\t", stringsAsFactors = F);
nsd_id <- match(node_short_df$pathway_name, all_nodes);
all_nodes[nsd_id] <- node_short_df$pathway_name_short;
# build data frame that contains node group info (how nodes are colored)
group_df <- read.delim(file = group_file, header = T, sep = "\t", stringsAsFactors = F);
all_node_group <- all_nodes;
names(all_node_group) <- all_nodes;
all_node_group[group_df$node_name] <- group_df$group;
node_df <- data.frame(node = all_nodes, group = all_node_group);
# use Sankey diagram to visualize CYP-related VNN connections
cyp_sk <- sankeyNetwork(
	Links = link_df,
	Nodes = node_df,
	Source = "source_id",
	Target = "target_id",
	Value = "value",
	NodeID = "node",
	NodeGroup = "group",
	fontSize = 32,
	fontFamily = 'Helvetica',
	nodeWidth = 30,
#	nodePadding = 40,
#	height = 700,
#	width = 2000,
	sinksRight = TRUE
)
cyp_sk

## 3. Visualize ABC transporter-related VNN connections in Sankey diagram  
# specify plotting parameters
abc_genes       <- c("ABCB1", "ABCG2");
path_rela_file  <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/abc_pathway_relations.tsv";
drug_atc_file   <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/tox21_10k_library_atc_code.tsv";
atc_manual_file <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/abc_drugs_atc_manual.tsv";
atc_class_file  <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/abc_drugs_atc_class.tsv"
group_file      <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/abc_node_group.tsv";
# extract ABC transporter-related VNN paths identified for HEK293 cytotoxicity outcome 
cd_id <- which(connection_df$gene_name %in% abc_genes);
ag_connection_df <- connection_df[cd_id, ];
# build data frame that contains ABC transporter-related connections between compounds and gene features 
abc_gene_ids <- lapply(abc_genes, function(ag){
	ag_id <- which(ag_connection_df$gene_name %in% ag);
	return(ag_id);
});
abc_gene_compounds <- lapply(abc_gene_ids, function(agi) unique(ag_connection_df$compound_name[agi]));
agc_len <- sapply(abc_gene_compounds, length);
gene_target <- unlist(mapply(function(x, y) rep(x, y), abc_genes, agc_len, SIMPLIFY = F));
compound_source <- unlist(abc_gene_compounds);
cg_df <- data.frame(compound_source, gene_target);
colnames(cg_df) <- c("compound_name", "target_name");
# build data frame that contains mapping bewteen compounds and their ATC class annotations   
drug_atc_df <- read.delim(file = drug_atc_file, header = T, sep = "\t", stringsAsFactors = F);
cd_atc_df <- merge(drug_atc_df, cg_df, by = "compound_name");
cd_atc3_df <- unique(cd_atc_df[, c("compound_name", "atc_id3")]);
atc_manual_df <- read.delim(file = atc_manual_file, header = T, sep = "\t", stringsAsFactors = F);
cad_id <- which(cd_atc3_df$compound_name %in% atc_manual_df$compound_name);
oad_id <- setdiff(1:nrow(cd_atc3_df), cad_id);
cd_atc_df1 <- rbind(cd_atc3_df[oad_id, ], atc_manual_df);
# build data frame that contains ABC transporter-related connections between ATC classes and gene features 
atc_class_df <- read.delim(file = atc_class_file, header = T, sep = "\t", stringsAsFactors = F);
cd_atc_class_df <- merge(cd_atc_df1, atc_class_df, by = "atc_id3");
cg_class_df <- merge(cg_df, cd_atc_class_df, by = "compound_name");
cg_class_df1 <- unique(cg_class_df[, c("atc_id3_class", "target_name")]);
colnames(cg_class_df1) <- c("source", "target");
cg_class_df1$value <- mapply(function(ccdaic, ccdtn){
	ccdaic_id <- which(cg_class_df$atc_id3_class %in% ccdaic);
	ccdtn_id <- which(cg_class_df$target_name %in% ccdtn);
	ccd_inter_len <- length(intersect(ccdaic_id, ccdtn_id));
	return(ccd_inter_len);
}, cg_class_df1$source, cg_class_df1$target);
# build data frame that contains ABC transporter-related connections between gene features and pathway modules 
abc_gene_pathways <- lapply(abc_gene_ids, function(agi) unique(ag_connection_df$pathway_name2[agi]));
agp_len <- sapply(abc_gene_pathways, length);
gene_source <- unlist(mapply(function(x, y) rep(x, y), abc_genes, agp_len, SIMPLIFY = F));
pathway_target <- unlist(abc_gene_pathways);
gp_df <- data.frame(gene_source, pathway_target);
colnames(gp_df) <- c("source", "target");
gp_df$value <- 1
abc_path_df <- rbind(cg_class_df1, gp_df);
# read in data frame that contains ABC transporter-related connections among pathway modules 
path_rela_df <- read.delim(file = path_rela_file, header = T, sep = "\t", stringsAsFactors = F);
path_rela_df$value <- 1
link_df <- rbind(abc_path_df, path_rela_df);
# build data frame that contains network node and edge info
all_nodes <- c(unique(cg_class_df1$source), abc_genes, unique(pathway_target), unique(path_rela_df$target));
link_df$source_id <- match(link_df$source, all_nodes) - 1;
link_df$target_id <- match(link_df$target, all_nodes) - 1;
# build data frame that contains node group info (how nodes are colored)
group_df <- read.delim(file = group_file, header = T, sep = "\t", stringsAsFactors = F);
all_node_group <- all_nodes;
names(all_node_group) <- all_nodes;
all_node_group[group_df$node_name] <- group_df$group;
node_df <- data.frame(node = all_nodes, group = all_node_group);
#  use Sankey diagram to visualize ABC transporter-related VNN connections
abc_sk <- sankeyNetwork(
	Links = link_df,
	Nodes = node_df,
	Source = "source_id",
	Target = "target_id",
	Value = "value",
	NodeID = "node",
	NodeGroup = "group",
	fontSize = 25,
	fontFamily = 'Helvetica',
	nodeWidth = 30,
#	nodePadding = 40,
#	height = 700,
#	width = 2000,
	sinksRight = TRUE
)
abc_sk

## 4. Visualize dasatinib-related VNN connections in Sankey diagram  
# specify plotting parameters
das_drug	<- "dasatinib";
path_rela_file  <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/dasatinib_pathway_relations.tsv";
path_short_file <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/dasatinib_pathway_short_names.tsv";
group_file      <- "data/compound_select_gene_comptoxai_tox21_interpret_analysis/dasatinib_node_group.tsv";
# extract dasatinib-related VNN paths identified for HEK293 cytotoxicity outcome 
cd_id <- which(connection_df$compound_name %in% das_drug);
cg_connection_df <- connection_df[cd_id, ];
# build data frame that contains dasatinib-related connections between compounds, gene features, and pathway modules 
cg_df <- unique(cg_connection_df[, c("compound_name", "gene_name")]);
gp_df <- unique(cg_connection_df[, c("gene_name", "pathway_name2")]);
colnames(cg_df) <- colnames(gp_df) <- c("source", "target");
das_path_df <- rbind(cg_df, gp_df);
# read in data frame that contains dasatinib-related connections among pathway modules
path_rela_df <- read.delim(file = path_rela_file, header = T, sep = "\t", stringsAsFactors = F);
link_df <- rbind(das_path_df, path_rela_df);
link_df$value <- 1;
# build data frame that contains network node and edge info
all_nodes <- c(das_drug, cg_df$target, unique(gp_df$target), unique(path_rela_df$target));
link_df$source_id <- match(link_df$source, all_nodes) - 1;
link_df$target_id <- match(link_df$target, all_nodes) - 1;
# build data frame that contains node name annotations
node_short_df <- read.delim(file = path_short_file, header = T, sep = "\t", stringsAsFactors = F);
nsd_id <- match(node_short_df$pathway_name, all_nodes);
all_nodes[nsd_id] <- node_short_df$pathway_name_short;
# build data frame that contains node group info (how nodes are colored)
group_df <- read.delim(file = group_file, header = T, sep = "\t", stringsAsFactors = F);
all_node_group <- all_nodes;
names(all_node_group) <- all_nodes;
all_node_group[group_df$node_name] <- group_df$group;
node_df <- data.frame(node = all_nodes, group = all_node_group);
# use Sankey diagram to visualize CYP-related VNN connections
das_sk <- sankeyNetwork(
	Links = link_df,
	Nodes = node_df,
	Source = "source_id",
	Target = "target_id",
	Value = "value",
	NodeID = "node",
	NodeGroup = "group",
	fontSize = 30,
	fontFamily = 'Helvetica',
	nodeWidth = 30,
	nodePadding = 40,
#	height = 700,
#	width = 2000,
	sinksRight = TRUE
)
das_sk


