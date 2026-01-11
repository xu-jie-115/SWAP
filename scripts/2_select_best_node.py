# ================================================================================
# Step 2: select best nodes (highest rank correlation)
# ================================================================================
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

min_data_points=10

df_swot = pd.read_csv('../data/swot_data.csv')

# preprocess SWOT data
node_counts = df_swot.groupby('node_id').size()
valid_nodes = node_counts[node_counts >= min_data_points].index
df_swot = df_swot[df_swot['node_id'].isin(valid_nodes)].copy()

dict_node = {}
for node_id, group in df_swot.groupby('node_id'):
    dict_node[node_id] = {
        'width': group['width'].values,
        'wse': group['wse'].values,
        'stationid': group['stationid'].iloc[0]
    }

node_ids = sorted(list(dict_node.keys()))
df_node = pd.DataFrame(index=node_ids, columns=['stationid', 'rank_corr'])
df_node.index.name = 'node_id'

# compute rank correlation (Spearman R)
for node_id in node_ids:
    node_data = dict_node.get(node_id)
    if node_data is None or len(node_data['width']) < 5:
        df_node.loc[node_id] = {
            'stationid': node_data['stationid'] if node_data else None,
            'rank_corr': 0.0
        }
        continue
    
    try:
        corr, _ = spearmanr(node_data['width'], node_data['wse'])
        if np.isnan(corr):
            corr = 0.0
    except:
        corr = 0.0
    
    df_node.loc[node_id] = {
        'stationid': node_data['stationid'],
        'rank_corr': corr
    }

df_node = df_node.dropna(subset=['stationid'])
df_node = df_node.reset_index()

# select best node for each station
df_node['rank_corr'] = pd.to_numeric(df_node['rank_corr'], errors='coerce')
max_idx = df_node.groupby('stationid')['rank_corr'].idxmax()
df_node_rmax = df_node.loc[max_idx]

# filter swot data
df_filtered = df_swot[df_swot['node_id'].isin(df_node_rmax['node_id'])]
df_filtered.to_csv('../results/2_swot_best_node.csv', index=False)