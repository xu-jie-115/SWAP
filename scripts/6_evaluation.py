# ================================================================================
# Step 6: evaluation (KGE, NSE, NRMSE)
# ================================================================================
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error

# NSE (Nash-Sutcliffe Efficiency)
def nse(observed, simulated):
    return 1 - (np.sum((observed - simulated) ** 2) / np.sum((observed - np.mean(observed)) ** 2))

# KGE (Kling-Gupta Efficiency)
def kge(observed, simulated):
    r = np.corrcoef(observed, simulated)[0, 1]
    alpha = np.mean(simulated) / np.mean(observed)
    beta  = np.std(simulated)/np.mean(simulated) / (np.std(observed)/np.mean(observed))
    return 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)

# Normalized RMSE
def normalized_rmse(observed, simulated):
    rmse = np.sqrt(mean_squared_error(observed, simulated))
    return rmse / np.mean(observed)

def fit_function(w,  z0, u1, s1):
    return z0 + u1 * (w ** s1)

df_fit_all = pd.read_csv('../results/4_fit_hypsometry.csv')
df_fit_all = df_fit_all[df_fit_all['fail'].isna()]  # use reliable results only
df_hyp_all = pd.read_csv('../results/5_hypsometry_median.csv')

stationids = sorted(df_fit_all['stationid'].unique())
df_res = []
for stationid in stationids:

    df_fit = df_fit_all[df_fit_all['stationid']==stationid]
    df_hyp = df_hyp_all[df_hyp_all['stationid']==stationid].reset_index(drop=True)

    # read discharge observation for evaluation
    df_eval = pd.read_csv('../data/discharge_obs/%s.csv'%stationid)
    df_eval = df_eval.rename(columns={'glow-mean': 'width'})
    df_eval.insert(0, 'stationid', stationid)
    w_low, w_high, slp = df_fit.iloc[0][['w_low','w_high','slp']]

    df_eval = df_eval[(df_eval['width']>=w_low) & (df_eval['width']<=w_high)]
    df_eval = df_eval.drop_duplicates('date')

    # interpolate XS area
    df_eval['idx_w'] = np.searchsorted(df_hyp['width'], df_eval['width'])
    df_eval['width_i-1'] = df_eval['idx_w'].apply(lambda x: df_hyp.loc[x-1,'width'])
    df_eval['width_i'] = df_eval['idx_w'].apply(lambda x: df_hyp.loc[x,'width'])
    df_eval['wse_i-1'] = df_eval['idx_w'].apply(lambda x: df_hyp.loc[x-1,'wse'])
    df_eval['wse_i'] = df_eval['idx_w'].apply(lambda x: df_hyp.loc[x,'wse'])
    df_eval['area_i-1'] = df_eval['idx_w'].apply(lambda x: df_hyp.loc[x-1,'area'])
    df_eval['area_hypso'] = df_eval.apply(lambda x: x['area_i-1'] + 0.5 * (x['width_i-1'] + x['width']) * \
        (x['wse_i'] - x['wse_i-1']) * (x['width'] - x['width_i-1']) / (x['width_i'] - x['width_i-1']), axis=1)
    
    # estimate discharge using Manning's equation
    df_eval['Q_est'] = df_eval['area_hypso']**(5.0/3.0) * df_eval['width']**(-2.0/3.0) * slp**0.5 / 0.035

    # compute KGE, NSE, NRMSE
    df_eval['kge'] = kge(df_eval['Q'], df_eval['Q_est'])
    df_eval['nse'] = nse(df_eval['Q'], df_eval['Q_est'])
    df_eval['nrmse'] = normalized_rmse(df_eval['Q'], df_eval['Q_est'])
    df_res.append(df_eval[['stationid','date','width','area_hypso','Q','Q_est','kge','nse','nrmse']])

df_res = pd.concat(df_res)
df_res.to_csv('../results/6_evaluation.csv', index=False)