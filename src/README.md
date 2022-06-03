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
  + [`predict_dtox.py`](predict_dtox.py) implements trained AIDTox model to predict outcome probability based on input feature data.
  + [`interpret_dtox.py`](`interpret_dtox.py`) implements layer-wise relevance propagation to evaluate relevance of AIDTox paths.
    + [`dtox_lrp.py`](dtox_lrp.py) contains functions used for implementing LRP to evaluate relevance of AIDTox paths.

+ Simple machine learning model
  + [`simple/simple.py`](simple/simple.py) develops and evaluates simple machine learning model (random forest or gradient boosting).
    + [`simple/simple_learning.py`](simple/simple_learning.py) contains functions for building, evaluating, and implementing simple machine learning models.
    + [`run/run_simple.R`](run/run_simple.R) generates shell scripts that run simple machine learning models on Tox21 datasets under different hyperparameter settings.

+ Model performance analysis, comparison, and visualization 
  + [`analysis_dtox/collect_model_results.R`](analysis_dtox/collect_model_results.R) collects machine learning model basic info and performance metrics from performance files.
  + [`analysis_dtox/analyze_dtox_results.py`](analysis_dtox/analyze_dtox_results.py) identifies optimal hyperparameter setting of machine learning method implementation, then compares and visualizes model performance across different method implementations.
    + [`analysis_dtox/dtox_analysis.py`](analysis_dtox/dtox_analysis.py) contains functions used in AIDTox model result anaysis.
    + [`analysis_dtox/dtox_plot.py`](analysis_dtox/dtox_plot.py) contains functions for visualizing AIDTox model results.

