# ================================================================================
# Step 5: compute median hypsometric curve for each station
# ================================================================================
import pandas as pd
import numpy as np

df_fit_all = pd.read_csv('../results/4_fit_hypsometry.csv')
df_fit_all = df_fit_all[df_fit_all['fail'].isna()]  # use reliable results only
stationids = sorted(df_fit_all['stationid'].unique())

df_res = []
for stationid in stationids:

    df_fit = df_fit_all[df_fit_all['stationid']==stationid].copy()
    w_low, w_high, w50, a50 = df_fit.loc[
        df_fit.index[0], 
        ['w_low','w_high','w50','a50']
    ]
    wse0 = df_fit['wse0'].values
    a = df_fit['a'].values
    b = df_fit['b'].values

    # approximate the median curve using 100 points
    w_list = np.linspace(w_low, w_high, 100)
    h_list, h_max, h_min = [], [], []  # median WSE, max WSE, min WSE among all curves
    for x in w_list:
        heights = wse0 + a * x**b
        h_list.append(np.median(heights))
        h_max.append(np.max(heights))
        h_min.append(np.min(heights))

    df_med = pd.DataFrame({
        'stationid': stationid,
        'width':w_list,
        'wse':h_list,
        'wse_max':h_max,
        'wse_min':h_min
    })
    
    # search index of median width (w50)
    idx50 = np.searchsorted(w_list, w50)
    # linear interpolation to derive WSE at median width
    h50 = (h_list[idx50-1] * (w_list[idx50]-w50) + h_list[idx50] * (w50-w_list[idx50-1])) / \
          (w_list[idx50]-w_list[idx50-1])
    
    # compute XS area based on the hypsometric curve
    # median XS area was derived in Step 4
    # areas higher than median are derived by adding Δarea to median XS area
    # areas lower than median are derived by substracting Δarea from median XS area
    df_med.loc[idx50, 'area'] = a50 + 0.5 * (w50 + w_list[idx50]) * (h_list[idx50] - h50)
    for i in np.arange(idx50+1, len(w_list)):
        df_med.loc[i,'area'] = df_med.loc[i-1,'area'] + 0.5 * (w_list[i-1]+w_list[i]) * (h_list[i]-h_list[i-1])
    for i in np.arange(idx50-1, -1, -1):
        df_med.loc[i,'area'] = df_med.loc[i+1,'area'] - 0.5 * (w_list[i+1]+w_list[i]) * (h_list[i+1]-h_list[i])
    df_res.append(df_med)

df_res = pd.concat(df_res)
df_res.to_csv('../results/5_hypsometry_median.csv', index=False)