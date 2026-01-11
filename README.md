# Leveraging "SWOT and A-Priori information (SWAP)" constrained channel parameters for improved historical river discharge estimates from space

## Overview
This repository provides a complete implementation of the methodology presented in our paper for improved river discharge estimation. The workflow processes Landsat-derived river width data and SWOT observations to generate hypsometric curves and evaluate discharge estimates against in-situ measurements.

## Repository Structure
```
├── scripts/
│   ├── 1_width_statistics.py          # Compute river width statistics from time series
│   ├── 2_select_best_node.py          # Select nodes with highest rank correlation
│   ├── 3_filter_swot.py               # Filter SWOT data (uncertainty, inversion, outlier)
│   ├── 4_fit_hypsometry.py            # Fit hypsometric curve with power function
│   ├── 5_derive_median_hypsometry.py  # Compute median hypsometric curves
│   └── 6_evaluation.py                # Evaluate against discharge observations
├── data/
│   ├── attributes.csv                 # River attributes
│   ├── glow_width_timeseries.csv      # Landsat-derived river width time series
│   ├── swot_data.csv                  # Preprocessed SWOT RiverSP Version D Node data
│   └── discharge_obs/                 # Discharge observations
├── results/
│   ├── 1_width_statistics.csv         # Width statistics results
│   ├── 2_swot_best_node.csv           # Selected nodes results
│   ├── 3_swot_filter.csv              # Filtered SWOT data
│   ├── 4_fit_hypsometry.csv           # Fitted hypsometric curves
│   ├── 5_hypsometry_median.csv        # Median hypsometric curves
│   └── 6_evaluation.csv               # Evaluation metrics (KGE, NSE, NRMSE)
└── README.md
```

## Script Descriptions
### 1. Width Statistics (`1_width_statistics.py`)
Computes river width statistics from Landsat-derived river width data, providing foundational width characteristics for subsequent analysis.

### 2. Node Selection (`2_select_best_node.py`)
Identifies optimal river nodes by selecting those with the highest rank correlation between width and WSE.

### 3. SWOT Filtering (`3_filter_swot.py`)
Applies quality control filters to SWOT data, removing records with high uncertainty, significant inversion, and outliers.

### 4. Hypsometry Fitting (`4_fit_hypsometry.py`)
Fits power-law hypsometric curves using various parameter combinations to establish width-WSE relationships.

### 5. Median Hypsometry (`5_derive_median_hypsometry.py`)
Computes median hypsometric curves for each station, providing robust width-WSE relationships.

### 6. Performance Evaluation (`6_evaluation.py`)
Evaluates derived discharge estimates against in-situ observations using Kling-Gupta Efficiency (KGE), Nash-Sutcliffe Efficiency (NSE), and Normalized Root Mean Square Error (NRMSE).

## Usage
Execute the scripts sequentially from 1 to 6. Each script reads necessary input files from the `data/` folder and writes results to the `results/` folder.

## Citation
Xu, J.#, Yuan, Z.#, Xu, Y., & Lin, P*. Leveraging "SWOT and A-Priori information (SWAP)" constrained channel parameters for improved historical river discharge estimates from space. (in revision)  
\# These authors contributed equally  
\* Corresponding Author: Peirong Lin (peironglinlin@pku.edu.cn)  
Contacts: Jie Xu (xujie115@pku.edu.cn), Zimin Yuan (ziminyuan@pku.edu.cn)