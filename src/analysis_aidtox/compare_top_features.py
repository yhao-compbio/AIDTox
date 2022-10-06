# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022
## This script visualizes the AIDTox model performance derived from different numbers of top predictive gene features, and compares the performance by connection types 


## Modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## 0. Input arguments 
perf_file	= 'data/compound_select_gene_comptoxai_tox21_implementation/compound_select_gene_comptoxai_tox21_implementation_optimal_performance_summary_by_training_root_loss.tsv'
class_file	= 'https://raw.githubusercontent.com/yhao-compbio/target/master/data/human_druggable_genome_protein_class.tsv'
feature_folder	= '/home/yunhao1/project/ontology/data/feature/comptoxai_select/comptoxai'
plot_folder	= 'plot/compound_select_gene_comptoxai_tox21_implementation/'
dataset_name	= ['tox21-rt-viability-hek293-p1', 'tox21-rt-viability-hepg2-p1']
dataset_title	= ['Cell viability (HEK293)', 'Cell viability (HepG2)']
connect_name	= ['binding', 'expression', 'combined']
connect_label	= ['binding', 'expression', 'hybrid']
feature_ns	= [100, 200, 300, 400, 500]

## 1. Use line charts to visualize the model performance ~ number of top features by connection type 
# read in the performance  
perf_df = pd.read_csv(perf_file, sep = '\t', header = 0)
# iterate by dataset 
perf_list = []
for ldn in range(0, len(dataset_name)):
	dn = dataset_name[ldn]
	# iterate by connection type
	dn_perf_list = []
	for cn in connect_name:
		# iterate by number of top features, obtain the index of performance metrics   
		cn_names = []	
		for fn in feature_ns:
			fn_name = cn + '_' + dn + '_' + str(fn)	
			cn_names.append(fn_name)
		# extract the performance metrics of models derived from current top features and current connection type
		cn_perf_df = perf_df[perf_df.dataset_name.isin(cn_names)]
		dn_perf_list.append(cn_perf_df)	
	perf_list.append(dn_perf_list)
	# specify figure and font size 
	fig = plt.figure(figsize = (6, 6))
	plt.rc('font', size = 20)
	plt.rc('axes', titlesize = 20)
	plt.rc('axes', labelsize = 20)
	plt.rc('xtick', labelsize = 20)
	plt.rc('ytick', labelsize = 20)
	plt.rc('legend', fontsize = 20)
	# iterate by connection type, plot model loss ~ number of top features 
	dn_perf_min = []
	dn_perf_max = []
	ax = fig.add_subplot(111)
	for ldpl in range(0, len(dn_perf_list)):
		# plot model loss ~ number of top features for the current connection type  
		dpf = dn_perf_list[ldpl]
		dpf_train_loss = dpf.training_root_loss.values
		ax.plot(feature_ns, dpf_train_loss, marker = 'o', label = connect_label[ldpl])
		# obtain the min and max value of model loss for the current connection type 
		dn_perf_min.append(np.min(dpf_train_loss))
		dn_perf_max.append(np.max(dpf_train_loss))
	# set x and y axis range 
	tl_min = np.min(dn_perf_min)
	tl_max = np.max(dn_perf_max)
	y_min = (np.ceil(tl_min * 50) - 1)/50
	y_max = (np.floor(tl_max * 50) + 1)/50	
	ax.set_ylim(y_min, y_max)
	# set x and y axis label
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlabel('# Top ranked gene features')
	ax.set_ylabel('Training BCE loss')
	# set figure title and legend 
	ax.set_title(dataset_title[ldn])
	ax.legend(frameon = False)
	# save figure
	fig.tight_layout()
	plt.savefig(plot_folder + 'feature_number_comparison/compound_select_gene_comptoxai_tox21_implementation_feature_number_training_loss_comparison_linechart_' + dn + '.pdf')
	plt.close()

## 2. Use barplot to visualize the optimal model performance by connection type  
# specify the optimal number of top features (400 in this case) 
select_id = 3
# specify plotting parameters  
compare_measures = ['testing_auc', 'testing_bac', 'testing_f1']
compare_measure_names = ['Area under ROC curve', 'Balanced accuracy', 'F1 score']
compare_measure_lims = [(0.4, 0.9), (0.4, 0.8), (0, 0.8)]
compare_color = ['dodgerblue', 'salmon', 'peachpuff']
# iterate by dataset  
for ldn in range(0, len(dataset_name)):
	# specify figure and font size 
	f = plt.figure(figsize = (2*len(compare_measures), 5))
	plt.rc('font', size = 18)
	plt.rc('axes', titlesize = 18)
	plt.rc('axes', labelsize = 18)
	plt.rc('xtick', labelsize = 16)
	plt.rc('ytick', labelsize = 15)
	plt.rc('legend', fontsize = 14)
	# iterate by performance metric of interest 
	for lcm in range(0, len(compare_measures)):
		# extract performance metric values and 95% CIs of different method implementations
		cm = compare_measures[lcm]
		cm_value = []
		cm_error = []
		for pll in perf_list[ldn]:
			cm_value.append(pll[cm].values[select_id])
			cm_error.append(pll[cm + '_ci'].values[select_id])
		x_pos = np.arange(0, len(cm_value))
		# make barplot showing the extracted performance metric values and 95% CIs 
		ax = f.add_subplot(int('1' + str(len(compare_measures)) + str(lcm+1)))
		ax.bar(x_pos, cm_value, yerr = cm_error, align = 'center', color = compare_color[lcm], ecolor = 'black', capsize = 2)
		# remove axes on the top, bottom, and right
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		# set y axis (performance metric) range and label
		ax.set_ylim(compare_measure_lims[lcm])
		ax.set_ylabel(compare_measure_names[lcm])
		# set x axis (different methods) ticks 
		ax.set_xticks(x_pos)
		ax.set_xticklabels(connect_label, rotation = 90)
	# add dataset name as title to the barplot   
	f.suptitle(dataset_title[ldn], size = 20)
	# save barplot 
	ldn_out = plot_folder + 'connect_type_comparison/compound_select_gene_comptoxai_tox21_implementation_connect_type_testing_performance_comparison_barplot_' + dataset_name[ldn] + '.pdf'
	plt.tight_layout()
	plt.savefig(ldn_out)
	plt.close()

