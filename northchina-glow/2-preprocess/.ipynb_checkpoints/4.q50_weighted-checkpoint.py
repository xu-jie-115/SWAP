import numpy as np
import pandas as pd
from copy import deepcopy
from datetime import datetime, timedelta
from netCDF4 import Dataset
import warnings
warnings.filterwarnings('ignore')

date_origin = datetime(1979, 1, 1)
df_date = pd.DataFrame()
df_date['date'] = pd.date_range(start='1979-01-01', end='2019-12-31')
df_date['month'] = df_date['date'].apply(lambda x: x.month)

df_comid = pd.read_csv('gauge3000_q50_slp.csv')
df_comid = df_comid.drop_duplicates(subset='stationid')
df_w_all = pd.read_csv('gauge3000_GLOW_width.csv')
df_w_all = df_w_all.merge(df_comid[['COMID','stationid']], on='stationid', how='left')
df_w_all = df_w_all.sort_values('COMID')
df_w_all['region'] = df_w_all['COMID'].apply(lambda x: int(str(x)[0]))
df_comid = df_comid.set_index('stationid')
#breakpoint()

df_res = pd.DataFrame(columns=['stationid','COMID','q50','q50_weighted'])
df_res['stationid'] = sorted(df_w_all['stationid'].unique())
df_res = df_res.set_index('stationid')
for region in range(1,9):
    print('\n... reading region %d/8'%(region), end='  ')
    df_w_r = df_w_all[df_w_all['region']==region]
    if len(df_w_r) == 0: continue
    f = Dataset('/shared1/RESEARCH_DATA/GRFR/output_pfaf_0%d_1979-2019.nc'%region)
    arr_comid = f.variables['rivid'][:]
    arr_qout = f.variables['Qout'][:]
    f.close()
    stations_r = sorted(df_w_r['stationid'].unique())

    for s, j in zip(stations_r, range(len(stations_r))):
        print('\r... processing %d/%d in region %d/8'%(j+1,len(stations_r),region), end='  ')

        # calculate weight from GLOW width distribution
        comid = df_comid.loc[s,'COMID']
        df_w = df_w_r[df_w_r['COMID']==comid]
        df_w['month'] = df_w['date'].apply(lambda x: int(x.split('-')[1]))
        df_w['weight'] = df_w.groupby('month')['month'].transform('count')
        df_weight = pd.DataFrame({'month': np.arange(1,13)})
        df_weight = df_weight.merge(df_w.drop_duplicates('month')[['month','weight']], on='month', how='left')
        df_weight = df_weight.dropna()

        # weighted q50
        idx_comid = np.where(arr_comid==comid)[0][0]
        qout = arr_qout[:,idx_comid]
        df_q = deepcopy(df_date)
        df_q['q'] = qout
        df_q = df_q.merge(df_weight, on='month', how='right')
        df_q = df_q.sort_values('q')
        df_q['cum_weight'] = df_q['weight'].cumsum()
        idx_median = np.where(df_q['cum_weight'].values<=df_q['weight'].sum()/2)[0][-1]
        q50_weighted = np.interp(x=df_q['weight'].sum()/2,
                                xp=df_q.iloc[idx_median:idx_median+2]['cum_weight'],
                                fp=df_q.iloc[idx_median:idx_median+2]['q'])
        df_res.loc[s,'COMID'] = comid
        df_res.loc[s,'q50'] = np.median(qout.data)
        df_res.loc[s,'q50_weighted'] = q50_weighted

df_res.to_csv('q50_weighted.csv')