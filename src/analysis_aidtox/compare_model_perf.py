# !/usr/bin/env python
## created by Yun Hao @MooreLab 2022
## This script visualizes the comparison of optimal performance between AIDTox and well-established models (DTox, RF, GB) using barplot.


## Modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


## 0. Input arguments 
model_loc_file	= 'data/compound_select_gene_comptoxai_tox21_implementation/compound_select_gene_comptoxai_tox21_compare_files.tsv'
plot_folder	= 'plot/compound_select_gene_comptoxai_tox21_implementation/'
dataset_name	= ['tox21-pxr-p1', 'tox21-rt-viability-hek293-p1', 'tox21-rt-viability-hepg2-p1']
dataset_title	= ['PXR agonist', 'Cell viability (HEK293)', 'Cell viability (HepG2)']

## 1. Extract the optimal performance of each model 
# obtain the file name index for each dataset   
full_name1 = []
full_name2 = []
for dn in dataset_name:
	full_name1.append('binding_' + dn + '_400')
	full_name2.append(dn + '_MultiSURF_400')
# iterate by model, extract the optimal performance of each model 
model_loc_df = pd.read_csv(model_loc_file, sep = '\t', header = 0)
model_perf_list = []
for mlds in range(0, model_loc_df.shape[0]):
	# read in the performance metrics of current model  
	perf_df = pd.read_csv(model_loc_df.perf_file[mlds], sep = '\t', header = 0)
	# for AIDTox model
	if model_loc_df.method[mlds] == 'AIDTox':
		mlds_perf_df = perf_df[perf_df.dataset_name.isin(full_name1)]
	# all the other models 
	else:
		mlds_perf_df = perf_df[perf_df.dataset_name.isin(full_name2)]
	model_perf_list.append(mlds_perf_df)

## 2. Use barplot to visualize the optimal performance by model 
# specify plotting parameters     
compare_measures = ['testing_auc', 'testing_bac', 'testing_f1']
compare_measure_names = ['Area under ROC curve', 'Balanced accuracy', 'F1 score']
compare_measure_lims = [(0.5, 0.9), (0.4, 0.8), (0, 0.8)]
compare_color = ['dodgerblue', 'salmon', 'peachpuff']
# iterate by dataset  
for ldn in range(0, len(dataset_name)):
	# specify figure and font size 
	f = plt.figure(figsize = (2*len(compare_measures), 5))
	plt.rc('font', size = 18)
	plt.rc('axes', titlesize = 18)
	plt.rc('axes', labelsize = 18)
	plt.rc('xtick', labelsize = 15)
	plt.rc('ytick', labelsize = 15)
	plt.rc('legend', fontsize = 14)
	# iterate by performance metric of interest 
	for lcm in range(0, len(compare_measures)):
		# extract performance metric values and 95% CIs of different models
		cm = compare_measures[lcm]
		cm_value = []
		cm_error = []
		for mpl in model_perf_list:
			cm_value.append(mpl[cm].values[ldn])
			cm_error.append(mpl[cm + '_ci'].values[ldn])
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
		ax.set_xticklabels(model_loc_df.method.values, rotation = 90)
	# add dataset name as title to the barplot   
	f.suptitle(dataset_title[ldn], size = 20)
	# save barplot 
	ldn_out = plot_folder + 'model_comparison/compound_select_gene_comptoxai_tox21_model_testing_performance_comparison_barplot_' + dataset_name[ldn] + '.pdf'
	plt.tight_layout()
	plt.savefig(ldn_out)
	plt.close()
