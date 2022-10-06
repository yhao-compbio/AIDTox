# This folder contains source code used by the repository.

## R/python scripts 

+ AIDTox model 
  + [`dtox.py`](dtox.py) learns and evaluates AIDTox model
    + [`dtox_data.py`](dtox_data.py) contains data-formatting functions used for AIDTox model training.
    + [`dtox_hierarchy.py`](dtox_hierarchy.py) contains functions used to process sorted AIDTox hiearchy files and compute model statistics.
    + [`dtox_nn.py`](dtox_nn.py) contains functions used to build basic neural network structure for AIDTox model.
    + [`dtox_loss.py`](dtox_loss.py) contains the the loss function used in AIDTox model.
    + [`early_stop.py`](early_stop.py) contains early stop function of AIDTox model.
    + [`dtox_learning.py`](dtox_learning.py) contains deep learning functions used in the AIDTox model construction.
    + [`run/run_aidtox_implementation.R`](run/run_aidtox_implementation.R) generates shell scripts that run AIDTox model (DTox with chemical-gene connections from ComptoxAI) on Tox21 datasets. [`run/run_dtox_comparison_implementation.R`](run/run_dtox_comparison_implementation.R) generates shell scripts that run DTox on Tox21 datasets, which matches the configuration of optimal AIDTox model, for performance comparison. [`run/run_aidtox_null_implementation.R`](run/run_aidtox_null_implementation.R) generates shell scripts that run AIDTox model on outcome-shuffled Tox21 datasets under Reactome pathway hierarchy.
  + [`predict_dtox.py`](predict_dtox.py) implements trained AIDTox model to predict outcome probability based on input feature data.
  + [`interpret_dtox.py`](`interpret_dtox.py`) implements layer-wise relevance propagation to evaluate relevance of AIDTox paths.
    + [`dtox_lrp.py`](dtox_lrp.py) contains functions used for implementing LRP to evaluate relevance of AIDTox paths.
  + [`run/run_interpret_aidtox.R`](run/run_interpret_aidtox.R) generates shell scripts that runs interpretation procedure on optimal AIDTox models trained for Tox21 datasets.

+ Simple machine learning model
  + [`simple/simple.py`](simple/simple.py) develops and evaluates simple machine learning model (random forest or gradient boosting).
    + [`simple/simple_learning.py`](simple/simple_learning.py) contains functions for building, evaluating, and implementing simple machine learning models.
    + [`run/run_simple.R`](run/run_simple.R) generates shell scripts that run simple machine learning models on Tox21 datasets under different hyperparameter settings.

+ Model performance analysis 
  + [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) collects machine learning model basic info and performance metrics from performance files.
  + [`analysis_dtox/analyze_dtox_results.py`](analysis_dtox/analyze_dtox_results.py) identifies optimal hyperparameter setting of machine learning method implementation, then compares and visualizes model performance across different method implementations.
    + [`analysis_dtox/dtox_analysis.py`](analysis_dtox/dtox_analysis.py) contains functions used in AIDTox model result anaysis.
    + [`analysis_dtox/dtox_plot.py`](analysis_dtox/dtox_plot.py) contains functions for visualizing AIDTox model results.

+ AIDTox feature analysis, model performance comparison, and visualization
  + [`analysis_aidtox/compare_top_features.py`](analysis_aidtox/compare_top_features.py) visualizes the AIDTox model performance derived from different numbers of top predictive gene features, and compares the performance by connection types.
  + [`analysis_aidtox/compare_model_perf.py`](analysis_aidtox/compare_model_perf.py) visualizes the comparison of optimal performance between AIDTox and well-established models (DTox, RF, GB) using barplot.
  + [`analysis_aidtox/visualize_protein_classes.py`](analysis_aidtox/visualize_protein_classes.py) uses Radar plots to visualize the class and subclass distribution comparison between AIDTox and DTox target features.
  + [`analysis_aidtox/analyze_compound_path_map.R`](analysis_aidtox/analyze_compound_path_map.R) extracts info of compounds and their identified VNN paths from AIDTox interpretation results.
  + [`analysis_aidtox/visualize_compound_path_map.R`](analysis_aidtox/visualize_compound_path_map.R) uses Sankey diagram to visualize VNN paths identified for HEK293 cytotoxicity outcome, which connect together compounds, gene features, pathway modules and the outcome.

+ [`functions.R`](functions.R) contains R functions required for other scripts in the repository.

## Executable shell scripts

+ AIDTox model implementation
  + [`run/run_aidtox_implementation_compound_gene_comptoxai_tox21.sh`](run/run_aidtox_implementation_compound_gene_comptoxai_tox21.sh) runs [`run/run_aidtox_implementation.R`](run/run_aidtox_implementation.R) to generate [`run/aidtox_compound_gene_comptoxai_tox21/dtox_select_compound_gene_comptoxai_tox21_implementation.sh`](run/aidtox_compound_gene_comptoxai_tox21/dtox_select_compound_gene_comptoxai_tox21_implementation.sh). [`run/aidtox_compound_gene_comptoxai_tox21/dtox_select_compound_gene_comptoxai_tox21_implementation.sh`](run/aidtox_compound_gene_comptoxai_tox21/dtox_select_compound_gene_comptoxai_tox21_implementation.sh) implements [`dtox.py`](dtox.py) on compound ComptoxAI connections-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy. 
  + [`run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null.sh`](run/aidtox_compound_gene_comptoxai_tox21_null/dtox_compound_select_gene_comptoxai_binding_tox21_null.sh) implements [`dtox.py`](dtox.py) on compound ComptoxAI connections-Tox21 assay outcome datasets under shuffled Reactome pathway hierarchy.

+ AIDTox model interpretation 
  + [`run/interpret_dtox_compound_select_gene_comptoxai_tox21_implementation.sh`](run/interpret_dtox_compound_select_gene_comptoxai_tox21_implementation.sh) implements [`interpret_dtox.py`](`interpret_dtox.py`) on optimal models trained for compound ComptoxAI connections-Tox21 assay outcome datasets.

+ Compared model implementation 
  + [`run/run_dtox_comparison_implementation_compound_gene_comptoxai_tox21.sh`](run/run_dtox_comparison_implementation_compound_gene_comptoxai_tox21.sh) runs [`run/run_dtox_comparison_implementation.R`](run/run_dtox_comparison_implementation.R) to generate [`run/dtox_comparison_compound_gene_comptoxai_tox21/dtox_compound_select_gene_comptoxai_target_probability_tox21_implementation.sh`](run/dtox_comparison_compound_gene_comptoxai_tox21/dtox_compound_select_gene_comptoxai_target_probability_tox21_implementation.sh). [`run/dtox_comparison_compound_gene_comptoxai_tox21/dtox_compound_select_gene_comptoxai_target_probability_tox21_implementation.sh`](run/dtox_comparison_compound_gene_comptoxai_tox21/dtox_compound_select_gene_comptoxai_target_probability_tox21_implementation.sh) implements [`dtox.py`](dtox.py) on compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy. 
  + [`run/run_simple_compound_structure_tox21_comptoxai.sh`](run/run_simple_compound_structure_tox21_comptoxai.sh) runs [`run/run_simple.R`](run/run_simple.R) to generate [`run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_randomforest.sh`](run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_randomforest.sh) and [`run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_xgboost.sh`](run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_xgboost.sh). [`run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_randomforest.sh`](run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_randomforest.sh) implements [`simple/simple.py`](simple/simple.py) to build random forest models on compound structure-Tox21 assay outcome datasets under different hyperparameter settings. [`run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_xgboost.sh`](run/simple_compound_structure_tox21/simple_compound_structure_tox21_select_gene_comptoxai_binding_xgboost.sh) implements [`simple/simple.py`](simple/simple.py) to build gradient boosting models on compound structure-Tox21 assay outcome datasets under different hyperparameter settings. 

+ Model performance analysis, comparison, and visualization
  + Result collection
    + [`run/collect_model_results_aidtox_compound_gene_comptoxai_tox21.sh`](run/collect_model_results_aidtox_compound_gene_comptoxai_tox21.sh) implements [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) to collect results of AIDTox models built upon compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.
    + [`run/collect_model_results_dtox_comparison_compound_gene_comptoxai_tox21.sh`](run/collect_model_results_dtox_comparison_compound_gene_comptoxai_tox21.sh) implements [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) to collect results of DTox comparison models built upon compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.
    + [`run/collect_model_results_simple_compound_structure_tox21.sh`](run/collect_model_results_simple_compound_structure_tox21.sh) implements [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) to collect results of simple machine learning models built upon compound structure-Tox21 assay outcome datasets under different hyperparameter settings.
  + Result analysis
    + [`run/analyze_dtox_results_aidtox_compound_gene_comptoxai_tox21.sh`](run/analyze_dtox_results_aidtox_compound_gene_comptoxai_tox21.sh) implements [`analysis_dtox/analyze_dtox_results.py`](analysis_dtox/analyze_dtox_results.py) to identify optimal hyperparameter setting of AIDTox implementation on compound target binding-Tox21 assay outcome datasets under sorted Reactome pathway hierarchy.
    + [`run/analyze_dtox_results_simple_compound_structure_tox21.sh`](run/analyze_dtox_results_simple_compound_structure_tox21.sh) implements [`analysis_dtox/analyze_dtox_results.py`](analysis_dtox/analyze_dtox_results.py) to identify optimal hyperparameter setting of simple machine learning model implementation on compound structure-Tox21 assay outcome datasets.
