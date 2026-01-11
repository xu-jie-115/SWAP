# ================================================================================
# Step 1: compute river width statistics from Landsat-derive river width
# ================================================================================
import pandas as pd
import numpy as np

min_width = 30

df_l8 = pd.read_csv('../data/glow_width_timeseries.csv')

stationids = df_l8['stationid'].unique()
result_data = []

for stationid in stationids:
    station_data_all = df_l8[df_l8['stationid'] == stationid]['width'].dropna()
    station_data = station_data_all[station_data_all >= min_width]  # require width >= 30m
    
    if len(station_data_all) == 0:
        continue
    if len(station_data) / len(station_data_all) < 0.95:
        continue
    
    if len(station_data) > 10:
        w50 = station_data.median()
        
        # compute Q1, Q3 and IQR
        q1 = station_data.quantile(0.25)
        q3 = station_data.quantile(0.75)
        iqr = q3 - q1
        
        # check whether width variablity is sufficiently large (IQR >= 5m)
        if iqr < 5:
            continue
        
        row_data = {
            'stationid': stationid,
            'w50': w50,
            'Q1': q1,
            'Q3': q3,
            'IQR': iqr
        }
        
        # compute river width at low/high flow
        w_low = max(q1 - 1.5 * iqr, min_width)  # w_low should be >= 30m
        w_high = q3 + 1.5 * iqr
        row_data['w_low'] = w_low
        row_data['w_high'] = w_high
        
        result_data.append(row_data)

df_w_stats = pd.DataFrame(result_data)
df_w_stats.to_csv('../results/1_width_statistics.csv', index=False)