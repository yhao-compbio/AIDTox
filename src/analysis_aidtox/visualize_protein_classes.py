# !/usr/bin/env python
## Created by Yun Hao @MooreLab 2022, adapted from https://matplotlib.org/stable/gallery/specialty_plots/radar_chart.html
## This script uses Radar plots to visualize the class and subclass distribution comparison between AIDTox and DTox target features  


## Modules
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D


## This function creates a RadarAxes projection and registers it
def radar_factory(num_vars, frame = 'circle'):
	## 0. Input arguments  
		# num_vars : int Number of variables for radar chart
		# frame : {'circle', 'polygon'}: Shape of frame surrounding axes.
	## 1. Calculate evenly-spaced axis angles
	theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)
	class RadarTransform(PolarAxes.PolarTransform):
		def transform_path_non_affine(self, path):
			# Paths with non-unit interpolation steps correspond to gridlines, in which case we force interpolation (to defeat PolarTransform's autoconversion to circular arcs).
			if path._interpolation_steps > 1:
				path = path.interpolated(num_vars)
			return Path(self.transform(path.vertices), path.codes)
	## 2. Define function classes of RadarAxes
	#
	class RadarAxes(PolarAxes):
		name = 'radar'
		# use 1 line segment to connect specified points
		RESOLUTION = 1
		PolarTransform = RadarTransform
		#
		def __init__(self, *args, **kwargs):
			super().__init__(*args, **kwargs)
			# rotate plot such that the first axis is at the top
			self.set_theta_zero_location('N')
		#
		def fill(self, *args, closed=True, **kwargs):
			"""Override fill so that line is closed by default"""
			return super().fill(closed=closed, *args, **kwargs)
		# 
		def plot(self, *args, **kwargs):
			"""Override plot so that line is closed by default"""
			lines = super().plot(*args, **kwargs)
			for line in lines:
				self._close_line(line)
		#
		def _close_line(self, line):
			x, y = line.get_data()
			# markers at x[0], y[0] get doubled-up
			if x[0] != x[-1]:
				x = np.append(x, x[0])
			y = np.append(y, y[0])
			line.set_data(x, y)
		#
		def set_varlabels(self, labels):
			self.set_thetagrids(np.degrees(theta), labels)
		#
		def _gen_axes_patch(self):
			# The Axes patch must be centered at (0.5, 0.5) and of radius 0.5 in axes coordinates.
			if frame == 'circle':
				return Circle((0.5, 0.5), 0.5)
			elif frame == 'polygon':
				return RegularPolygon((0.5, 0.5), num_vars, radius=.5, edgecolor="k")
			else:
				raise ValueError("Unknown value for 'frame': %s" % frame)
		# 
		def _gen_axes_spines(self):
			if frame == 'circle':
				return super()._gen_axes_spines()
			elif frame == 'polygon':
				# spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
				spine = Spine(axes=self, spine_type='circle', path=Path.unit_regular_polygon(num_vars))
				# unit_regular_polygon gives a polygon of radius 1 centered at (0, 0) but we want a polygon of radius 0.5 centered at (0.5, 0.5) in axes coordinates.
				spine.set_transform(Affine2D().scale(.5).translate(.5, .5) + self.transAxes)
				return {'polar': spine}
			else:
				raise ValueError("Unknown value for 'frame': %s" % frame)
	register_projection(RadarAxes)
	return theta


## 0. Input arguments 
class_file	= '/home/yunhao1/project/target/data/human_druggable_genome_protein_class_v2.tsv'
feature_file1	= '/home/yunhao1/project/ontology/data/feature/target_profile_from_maccs_fingerprints.txt' 	
feature_file2	= '/home/yunhao1/project/ontology/data/feature/comptoxai_binding_genes.txt'
plot_folder     = 'plot/compound_select_gene_comptoxai_tox21_implementation/binding_protein_class/'

## 1. Obtain the protein class annotation of DTox and AIDTox input features 
# read in the class annotation of proteins in human druggable genome
protein_class = pd.read_csv(class_file, sep = '\t', header = 0)
protein_class = protein_class[protein_class['class'] != 'Other protein']
# read in input target features of DTox 
feature_df1 = pd.read_csv(feature_file1, sep = '\t', header = None)
feature_df1.columns = ['target']
# obtaion the class annotation of DTox target features
target_group_df1 = pd.merge(feature_df1, protein_class, how = 'left', left_on = 'target', right_on = 'uniprot_id')
target_group_df1 = target_group_df1.fillna('Others')
# read in input target gene features of AIDTox
feature_df2 = pd.read_csv(feature_file2, sep = '\t', header = None)
feature_df2.columns = ['target']
# obtaion the class annotation of AIDTox target features
target_group_df2 = pd.merge(feature_df2, protein_class, how = 'left', left_on = 'target', right_on = 'uniprot_id')
target_group_df2 = target_group_df2.fillna('Others')

## 2. Visualize the class distribution comparison between AIDTox and DTox target features by Radar plot 
# specify the target classes to be plotted  
class_od = ['Enzyme', 'G protein coupled receptor', 'Ion channel', 'Others', 'Transporter', 'Nuclear hormone receptor', 'Catalytic receptor']
class_label = ['Enzyme', 'G protein\n coupled\n receptor', 'Ion\n channel', 'Others', 'Transporter', 'Nuclear\n hormone\n receptor', 'Catalytic\n receptor']
# compute class distribution of DTox target features
target_group_count1 = target_group_df1.groupby('class')['target'].nunique()
count_df1 = target_group_count1.loc[class_od,]
# compute class distribution of AIDTox target features
target_group_count2 = target_group_df2.groupby('class')['target'].nunique()
count_df2 = target_group_count2.loc[class_od,]
# compute hypothetical distribution of AIDTox target features give a proportional increase from DTox 
count_df3 = count_df1 * count_df2.sum()/count_df1.sum()
# create a RadarAxes projection for the Radar plot 
N = len(class_od)
theta = radar_factory(N, frame = 'polygon')
# specify figure and font size 
plt.rc('font', size = 20)
plt.rc('axes', titlesize = 20)
plt.rc('axes', labelsize = 20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)
fig, ax = plt.subplots(figsize=(6, 6), nrows = 1, ncols = 1, subplot_kw = dict(projection='radar'))
fig.subplots_adjust(wspace = 0.25, hspace = 0.20, top = 0.8, bottom = 0.05)
# use radar plot to visualize the class distribution comparison between AIDTox and DTox target features 
ax.set_rgrids([100, 200, 300])
ax.plot(theta, count_df1.values, color = 'b')
ax.fill(theta, count_df1.values, facecolor = 'b', alpha = 0.25, label = '_nolegend_')
ax.plot(theta, count_df2.values, color = 'r')
ax.fill(theta, count_df2.values, facecolor = 'r', alpha = 0.25, label = '_nolegend_')
ax.plot(theta, count_df3.values, '--', color = 'tab:gray')
ax.set_varlabels(class_label)
# add figure legend 
labels = ('DTox', 'AIDTox')
legend = ax.legend(labels, loc = (0.7, 0.98), labelspacing = 0.1)
# save plot
plt.tight_layout()
plt.savefig(plot_folder + 'DTox_AIDTox_feature_space_protein_class_compare_radarchart.pdf')
plt.close()

## 3. Visualize the enzyme subclass distribution comparison between AIDTox and DTox target features by Radar plot
# specify the enzyme subclasses to be plotted  
subclass_od = ['Protease', 'PTP', 'CYP', 'PI3K', 'PI']
subclass_label = ['Protease', 'Protein\n tyrosine\n phosphatase', 'Cytochrome\n P450 enzyme', 'Phosphoinositide\n 3-kinase', 'Protease\n inhibitor']
# compute enzyme subclass distribution of DTox target features
target_group_count4 = target_group_df1.groupby('subclass')['target'].nunique()
count_df4 = target_group_count4.loc[subclass_od,]
# compute enzyme subclass distribution of AIDTox target features
target_group_count5 = target_group_df2.groupby('subclass')['target'].nunique()
count_df5 = target_group_count5.loc[subclass_od,]
# create a RadarAxes projection for the Radar plot 
N = len(subclass_od)
theta = radar_factory(N, frame = 'polygon')
# specify figure and font size
plt.rc('font', size = 20)
plt.rc('axes', titlesize = 20)
plt.rc('axes', labelsize = 20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('legend', fontsize = 20)
fig, ax = plt.subplots(figsize=(6, 6), nrows = 1, ncols = 1, subplot_kw = dict(projection='radar'))
fig.subplots_adjust(wspace = 0.25, hspace = 0.20, top = 0.8, bottom = 0.05)
# use radar plot to visualize the enzyme subclass distribution comparison between AIDTox and DTox target features 
ax.set_rgrids([20, 40])
ax.plot(theta, count_df4.values, color = 'b')
ax.fill(theta, count_df4.values, facecolor = 'b', alpha = 0.25, label = '_nolegend_')
ax.plot(theta, count_df5.values, color = 'r')
ax.fill(theta, count_df5.values, facecolor = 'r', alpha = 0.25, label = '_nolegend_')
ax.set_varlabels(subclass_label)
# add figure legend 
labels = ('DTox', 'AIDTox')
legend = ax.legend(labels, loc = (0.7, 0.98), labelspacing = 0.1)
# save plot
plt.tight_layout()
plt.savefig(plot_folder + 'DTox_AIDTox_feature_space_enzyme_subclass_compare_radarchart.pdf')
plt.close()
