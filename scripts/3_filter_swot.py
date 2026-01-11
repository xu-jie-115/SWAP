# ================================================================================
# Step 3: filter SWOT data (uncertainty, inversion, outlier)
# ================================================================================
import pandas as pd
import numpy as np

inverse_ratio_thresh = 0.5

df_w_stats = pd.read_csv('../results/1_width_statistics.csv').set_index('stationid')
df_swot = pd.read_csv('../results/2_swot_best_node.csv')

stationids = sorted(df_swot['stationid'].unique())
df_res = []

for stationid in stationids:
    if stationid not in df_w_stats.index:
        continue

    df = df_swot[df_swot['stationid'] == stationid].copy()
    if len(df) < 5:
        continue
    
    # Filter 1: uncertainty
    df['width_u_r'] = df['width_u'] / df['width']  # relative width uncertainty
    df_f1 = df[(df['wse_u'] <= 0.4) & (df['width_u_r'] <= 0.1)]
    if len(df_f1) < 5:
        continue
    
    # Filter 2: inversion
    df_f2 = df_f1.copy()
    while True:
        n = len(df_f2)
        w = df_f2['width'].values
        h = df_f2['wse'].values
        w_diff = np.repeat(w,n).reshape(n,-1) - np.tile(w,(n,1))
        h_diff = np.repeat(h,n).reshape(n,-1) - np.tile(h,(n,1))
        wh = w_diff * h_diff  # inversion matrix, value < 0 indicates inversion
        inverse = np.count_nonzero(wh<0, axis=1)
        idx_max = np.argmax(inverse)

        # stop filtering when proportion of inversion < 50% for all points
        if inverse[idx_max] / n < inverse_ratio_thresh: break
        df_f2 = df_f2.drop(df_f2.index[idx_max])

    if len(df_f2) < 5:
        continue
    
    # Filter 3: outlier
    w_low = df_w_stats.loc[stationid, 'w_low']
    w_high = df_w_stats.loc[stationid, 'w_high']
    d_bankfull = 0.27 * (w_high / 7.2) ** 0.6  # bankfull depth estimated from Andreadis (2013)
    h50 = df_f2['wse'].median()
    
    df_f3 = df_f2.copy()
    df_f3 = df_f3[(df_f3['width'] <= w_high) & (df_f3['width'] >= w_low)]  # remove width outliers
    df_f3 = df_f3[(df_f3['wse'] <= h50 + d_bankfull) & (df_f3['wse'] >= h50 - d_bankfull)]  # remove WSE outliers
    
    if len(df_f3) < 5:
        continue

    df_res.append(df_f3)

df_res = pd.concat(df_res)
df_res = df_res.drop_duplicates(subset=['node_id', 'date', 'stationid'])
df_res.to_csv('../results/3_swot_filter.csv', index=False)