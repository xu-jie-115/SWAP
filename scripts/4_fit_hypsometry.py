# ================================================================================
# Step 4: fit hypsometric curve with power function using diffrent parameters
# ================================================================================
import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.optimize import least_squares

# wse = wse0 + a * width**b
# return array of redisuals
def PowerFunc(params, X, y):
    wse0, a, b = params
    return y - (wse0 + a * X**b)

# soft_l1, z=residual^2
# 'weight' for swot points, low and high point equally divide remaining weight
def Loss(z):
    rho = np.zeros(3*len(z)).reshape(3,-1)
    rho[0] = 2 * ((1 + z)**0.5 - 1)
    rho[1] = (1 + z)**(-0.5)
    rho[2] = -0.5 * (1 + z)**(-1.5)
    rho[:,0] *= (len(z) - 2) / weight * (1 - weight) / 2  # (1-SWOT_weight)/2 for low point
    rho[:,1] *= (len(z) - 2) / weight * (1 - weight) / 2  # (1-SWOT_weight)/2 for high point
    return rho


R_list = np.array([0.5, 1, 2, 4, 8])
GAP_list = np.array([-0.3, -0.1, 0, 0.1, 0.3])
W_list = np.array([0.3, 0.5, 0.7])

df_swot = pd.read_csv('../results/3_swot_filter.csv')
df_attrs = pd.read_csv('../data/attributes.csv').set_index('COMID')
df_w_stats = pd.read_csv('../results/1_width_statistics.csv').set_index('stationid')

# prepare fitting result dataframe
df_res = df_swot.drop_duplicates('stationid')
df_res = df_res[['stationid','COMID']]
df_res = pd.concat([df_res] * (len(R_list) * len(GAP_list) * len(W_list)), ignore_index=True)
df_res = df_res.sort_values('stationid')
stationids = df_res['stationid'].unique()
df_res['R'] = R_list[np.tile(np.repeat(np.arange(len(R_list)), len(GAP_list)*len(W_list)), len(stationids))]
df_res['GAP'] = GAP_list[np.tile(np.repeat(np.arange(len(GAP_list)), len(W_list)), len(R_list)*len(stationids))]
df_res['W'] = W_list[np.tile(np.arange(len(W_list)), len(GAP_list)*len(R_list)*len(stationids))]
df_res = df_res.set_index(['stationid','R','GAP','W'])

q50, w50, slp, weight = None, None, None, None
for stationid in stationids[:]:

    print('... processing %s ...'%stationid)

    # extract data from files
    df = df_swot[df_swot['stationid']==stationid].copy()
    comid = df.loc[df.index[0], 'COMID']
    q50 = df_attrs.loc[comid, 'q50_weighted']  # median discharge
    slp = df_attrs.loc[comid, 'slope']  # channel slope
    w50, w_low, w_high = df_w_stats.loc[stationid, ['w50','w_low','w_high']]  # median width

    d_bankfull = 0.27 * (w_high / 7.2)**0.6
    swot_wsemax = df.sort_values('wse', ascending=False).iloc[0]
    d_wsemax = 0.27 * (swot_wsemax['width'] / 7.2)**0.6

    # compute median WSE (h50)
    # linear regression using 5 nearest SWOT data around median width (w50)
    df['w50_diff'] = np.abs(df['width'] - w50)
    df = df.sort_values('w50_diff')
    xdata = df.iloc[:5]['width'].values
    ydata = df.iloc[:5]['wse'].values
    res = linregress(xdata, ydata)
    if res[0] >= 0:
        h50 = res[0] * w50 + res[1]
    else:
        h50 = df.iloc[:5]['wse'].mean()

    # compute median XS area (a50) based on Manning's equation
    # solve (a50 / w50)**(2/3) * slp**(1/2) / 0.035 * a50 = q50
    a50 = (q50 * 0.035 / slp**(1/2) * w50**(2/3))**(3/5)

    # fit hypsometric curves with combinations of parameters
    for r_low in R_list:
        for gap in GAP_list:
            for weight in W_list:
                # compute WSE at low/high flow
                # suppose channel below median can be approximated by d = a * w**R
                # solve a50 = a * w50**R * w50 - a / (R+1) * w50**(R+1)
                a_low = a50 * (r_low+1) / r_low / w50**(r_low+1)
                h0 = h50 - a_low * w50**r_low
                h_low = h0 + a_low * w_low**r_low
                h_high = df['wse'].max() + (d_bankfull - d_wsemax) + gap * d_bankfull

                # fit
                xdata = df['width'].values
                ydata = df['wse'].values
                a_default = (h_high - h0) / w_high**2
                xdata = np.insert(xdata, 0, [w_low,w_high])
                ydata = np.insert(ydata, 0, [h_low,h_high])

                ls = least_squares(
                    PowerFunc,
                    x0=[h0,a_default,2],
                    loss=Loss,
                    args=(xdata,ydata),
                    max_nfev=100
                )

                df_res.loc[(stationid,r_low,gap,weight), ['wse0','a','b']] = ls.x
                df_res.loc[(stationid,r_low,gap,weight), ['a50','w50','q50']] = [a50, w50, q50]
                df_res.loc[(stationid,r_low,gap,weight), ['w_low','w_high','h_low','h_high','slp']] = \
                    [w_low, w_high, h_low, h_high, slp]

                # check if fitting result is reliable
                if ls.status <= 0:
                    # fitting failed (max iteration exceeded)
                    df_res.loc[(stationid,r_low,gap,weight), 'fail'] = 'ls_failed'
                elif ls.x[1] * ls.x[2] < 0:
                    # abnormal fitting results (a * b < 0)
                    df_res.loc[(stationid,r_low,gap,weight), 'fail'] = 'ab_negative'

df_res.to_csv('../results/4_fit_hypsometry.csv')
print('done')