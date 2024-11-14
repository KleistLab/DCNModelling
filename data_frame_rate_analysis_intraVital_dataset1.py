import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.style.use('default')


#Inverse logistic function
def func(x, L ,x0, k, b):
    y = 1/(L / (1 + np.exp(-k*(x-x0)))+b)
    return y

# Import data
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'DATA-set_Analyse_Intravital.xlsx')
data = pd.read_excel(filename, sheet_name = 'RAW_20170817_Ctrl')

# Keep only multi-fibre structures
data_multi = data[data.MULTI_TRACKID == 'MULTI']

# Keep only single-fibre structures
data_single = data[data.MULTI_TRACKID == 'SINGLE']

data_multi_unique  = data_multi.drop_duplicates('TrackID')

data_single_unique  = data_single.drop_duplicates('TrackID')

# Convert hours to minutes
data_multi["Time"] = data_multi["Time (h)"] * 60
data_single["Time"] = data_single["Time (h)"] * 60


# Total number of axonal structuress that had more than one fibre (multifibre) at some point in time
total_nr_multifibre_bundles = len(data_multi_unique)
list_trackID_multifibres_baseline = data_multi_unique['TrackID']

perc_total_nr_multifibre_bundles = total_nr_multifibre_bundles / total_nr_multifibre_bundles

# Total number of axonal structuress that had more only one fibre (single-fibre) at all time points
total_nr_singlefibre_bundles = len(data_single_unique)
list_trackID_singlefibres_baseline = data_single_unique['TrackID']

perc_total_nr_singlefibre_bundles = total_nr_singlefibre_bundles / total_nr_singlefibre_bundles


# %% Multifibre axonal structures

# Times are now from 0 to 1290 minutes, with a 30 minute frame rate


########### Considering permutations #######################

# data we get if we measure every 60 minutes

#permutation 1
times_60_perm1 = np.arange(0,1291,60)
data_multi_60min_rate = data_multi[data_multi['Time'].isin(times_60_perm1)] 
unique_multifibres_60min_rate  = data_multi_60min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_60min = unique_multifibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_multifibres_60min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_60min_rate_perm1 = len(IDs60min_inBaseline)
perc_multifibre_bundles_60min_rate_perm1 = nr_multifibre_bundles_60min_rate_perm1 / total_nr_multifibre_bundles

#permutation 2
times_60_perm2 = np.arange(30,1291,60)
data_multi_60min_rate = data_multi[data_multi['Time'].isin(times_60_perm2)] 
unique_multifibres_60min_rate  = data_multi_60min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_60min = unique_multifibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_multifibres_60min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_60min_rate_perm2 = len(IDs60min_inBaseline)
perc_multifibre_bundles_60min_rate_perm2 = nr_multifibre_bundles_60min_rate_perm2 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_60min_rate = np.mean([nr_multifibre_bundles_60min_rate_perm1, nr_multifibre_bundles_60min_rate_perm2])
perc_multifibre_bundles_60min_rate = np.mean([perc_multifibre_bundles_60min_rate_perm1, perc_multifibre_bundles_60min_rate_perm2])



# data we get if we measure every 90 minutes

#permutation 1
times_90_perm1 = np.arange(0,1291,90)
data_multi_90min_rate = data_multi[data_multi['Time'].isin(times_90_perm1)] 
unique_multifibres_90min_rate  = data_multi_90min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_90min = unique_multifibres_90min_rate['TrackID']
IDs90min_inBaseline = list_trackID_multifibres_90min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_90min_rate_perm1 = len(IDs90min_inBaseline)
perc_multifibre_bundles_90min_rate_perm1 = nr_multifibre_bundles_90min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_90_perm2 = np.arange(30,1291,90)
data_multi_90min_rate = data_multi[data_multi['Time'].isin(times_90_perm2)] 
unique_multifibres_90min_rate  = data_multi_90min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_90min = unique_multifibres_90min_rate['TrackID']
IDs90min_inBaseline = list_trackID_multifibres_90min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_90min_rate_perm2 = len(IDs90min_inBaseline)
perc_multifibre_bundles_90min_rate_perm2 = nr_multifibre_bundles_90min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_90_perm3 = np.arange(60,1291,90)
data_multi_90min_rate = data_multi[data_multi['Time'].isin(times_90_perm3)] 
unique_multifibres_90min_rate  = data_multi_90min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_90min = unique_multifibres_90min_rate['TrackID']
IDs90min_inBaseline = list_trackID_multifibres_90min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_90min_rate_perm3 = len(IDs90min_inBaseline)
perc_multifibre_bundles_90min_rate_perm3 = nr_multifibre_bundles_90min_rate_perm3 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_90min_rate = np.mean([nr_multifibre_bundles_90min_rate_perm1, nr_multifibre_bundles_90min_rate_perm2, 
                                            nr_multifibre_bundles_90min_rate_perm3])
perc_multifibre_bundles_90min_rate = np.mean([perc_multifibre_bundles_90min_rate_perm1, perc_multifibre_bundles_90min_rate_perm2, 
                                              perc_multifibre_bundles_90min_rate_perm3])



# data we get if we measure every 120 minutes

#permutation 1
times_120_perm1 = np.arange(0,1291,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm1)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm1 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm1 = nr_multifibre_bundles_120min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_120_perm2 = np.arange(30,1291,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm2)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm2 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm2 = nr_multifibre_bundles_120min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_120_perm3 = np.arange(60,1291,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm3)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm3 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm3 = nr_multifibre_bundles_120min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_120_perm4 = np.arange(90,1291,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm4)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm4 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm4 = nr_multifibre_bundles_120min_rate_perm4 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_120min_rate = np.mean([nr_multifibre_bundles_120min_rate_perm1, nr_multifibre_bundles_120min_rate_perm2, 
                                            nr_multifibre_bundles_120min_rate_perm3, nr_multifibre_bundles_120min_rate_perm4])
perc_multifibre_bundles_120min_rate = np.mean([perc_multifibre_bundles_120min_rate_perm1, perc_multifibre_bundles_120min_rate_perm2, 
                                              perc_multifibre_bundles_120min_rate_perm3, perc_multifibre_bundles_120min_rate_perm4])


# data we get if we measure every 150 minutes

#permutation 1
times_150_perm1 = np.arange(0,1291,150)
data_multi_150min_rate = data_multi[data_multi['Time'].isin(times_150_perm1)] 
unique_multifibres_150min_rate  = data_multi_150min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_150min = unique_multifibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_multifibres_150min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_150min_rate_perm1 = len(IDs150min_inBaseline)
perc_multifibre_bundles_150min_rate_perm1 = nr_multifibre_bundles_150min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_150_perm2 = np.arange(30,1291,150)
data_multi_150min_rate = data_multi[data_multi['Time'].isin(times_150_perm2)] 
unique_multifibres_150min_rate  = data_multi_150min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_150min = unique_multifibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_multifibres_150min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_150min_rate_perm2 = len(IDs150min_inBaseline)
perc_multifibre_bundles_150min_rate_perm2 = nr_multifibre_bundles_150min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_150_perm3 = np.arange(60,1291,150)
data_multi_150min_rate = data_multi[data_multi['Time'].isin(times_150_perm3)] 
unique_multifibres_150min_rate  = data_multi_150min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_150min = unique_multifibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_multifibres_150min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_150min_rate_perm3 = len(IDs150min_inBaseline)
perc_multifibre_bundles_150min_rate_perm3 = nr_multifibre_bundles_150min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_150_perm4 = np.arange(90,1291,150)
data_multi_150min_rate = data_multi[data_multi['Time'].isin(times_150_perm4)] 
unique_multifibres_150min_rate  = data_multi_150min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_150min = unique_multifibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_multifibres_150min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_150min_rate_perm4 = len(IDs150min_inBaseline)
perc_multifibre_bundles_150min_rate_perm4 = nr_multifibre_bundles_150min_rate_perm4 / total_nr_multifibre_bundles

#permuation 5
times_150_perm5 = np.arange(120,1291,150)
data_multi_150min_rate = data_multi[data_multi['Time'].isin(times_150_perm5)] 
unique_multifibres_150min_rate  = data_multi_150min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_150min = unique_multifibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_multifibres_150min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_150min_rate_perm5 = len(IDs150min_inBaseline)
perc_multifibre_bundles_150min_rate_perm5 = nr_multifibre_bundles_150min_rate_perm5 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_150min_rate = np.mean([nr_multifibre_bundles_150min_rate_perm1, nr_multifibre_bundles_150min_rate_perm2, 
                                            nr_multifibre_bundles_150min_rate_perm3, nr_multifibre_bundles_150min_rate_perm4,
                                            nr_multifibre_bundles_150min_rate_perm5])
perc_multifibre_bundles_150min_rate = np.mean([perc_multifibre_bundles_150min_rate_perm1, perc_multifibre_bundles_150min_rate_perm2, 
                                              perc_multifibre_bundles_150min_rate_perm3, perc_multifibre_bundles_150min_rate_perm4,
                                              perc_multifibre_bundles_150min_rate_perm5])

# data we get if we measure every 180 minutes

#permutation 1
times_180_perm1 = np.arange(0,1291,180)
data_multi_180min_rate = data_multi[data_multi['Time'].isin(times_180_perm1)] 
unique_multifibres_180min_rate  = data_multi_180min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_180min = unique_multifibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_multifibres_180min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_180min_rate_perm1 = len(IDs180min_inBaseline)
perc_multifibre_bundles_180min_rate_perm1 = nr_multifibre_bundles_180min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_180_perm2 = np.arange(30,1291,180)
data_multi_180min_rate = data_multi[data_multi['Time'].isin(times_180_perm2)] 
unique_multifibres_180min_rate  = data_multi_180min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_180min = unique_multifibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_multifibres_180min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_180min_rate_perm2 = len(IDs180min_inBaseline)
perc_multifibre_bundles_180min_rate_perm2 = nr_multifibre_bundles_180min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_180_perm3 = np.arange(60,1291,180)
data_multi_180min_rate = data_multi[data_multi['Time'].isin(times_180_perm3)] 
unique_multifibres_180min_rate  = data_multi_180min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_180min = unique_multifibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_multifibres_180min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_180min_rate_perm3 = len(IDs180min_inBaseline)
perc_multifibre_bundles_180min_rate_perm3 = nr_multifibre_bundles_180min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_180_perm4 = np.arange(90,1291,180)
data_multi_180min_rate = data_multi[data_multi['Time'].isin(times_180_perm4)] 
unique_multifibres_180min_rate  = data_multi_180min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_180min = unique_multifibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_multifibres_180min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_180min_rate_perm4 = len(IDs180min_inBaseline)
perc_multifibre_bundles_180min_rate_perm4 = nr_multifibre_bundles_180min_rate_perm4 / total_nr_multifibre_bundles

#permuation 5
times_180_perm5 = np.arange(120,1291,180)
data_multi_180min_rate = data_multi[data_multi['Time'].isin(times_180_perm5)] 
unique_multifibres_180min_rate  = data_multi_180min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_180min = unique_multifibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_multifibres_180min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_180min_rate_perm5 = len(IDs180min_inBaseline)
perc_multifibre_bundles_180min_rate_perm5 = nr_multifibre_bundles_180min_rate_perm5 / total_nr_multifibre_bundles

#permuation 6
times_180_perm6 = np.arange(150,1291,180)
data_multi_180min_rate = data_multi[data_multi['Time'].isin(times_180_perm6)] 
unique_multifibres_180min_rate  = data_multi_180min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_180min = unique_multifibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_multifibres_180min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_180min_rate_perm6 = len(IDs180min_inBaseline)
perc_multifibre_bundles_180min_rate_perm6 = nr_multifibre_bundles_180min_rate_perm6 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_180min_rate = np.mean([nr_multifibre_bundles_180min_rate_perm1, nr_multifibre_bundles_180min_rate_perm2, 
                                            nr_multifibre_bundles_180min_rate_perm3, nr_multifibre_bundles_180min_rate_perm4,
                                            nr_multifibre_bundles_180min_rate_perm5, nr_multifibre_bundles_180min_rate_perm6])
perc_multifibre_bundles_180min_rate = np.mean([perc_multifibre_bundles_180min_rate_perm1, perc_multifibre_bundles_180min_rate_perm2, 
                                              perc_multifibre_bundles_180min_rate_perm3, perc_multifibre_bundles_180min_rate_perm4,
                                              perc_multifibre_bundles_180min_rate_perm5, perc_multifibre_bundles_180min_rate_perm6])


# data we get if we measure every 210 minutes

#permutation 1
times_210_perm1 = np.arange(0,1291,210)
data_multi_210min_rate = data_multi[data_multi['Time'].isin(times_210_perm1)] 
unique_multifibres_210min_rate  = data_multi_210min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_210min = unique_multifibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_multifibres_210min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_210min_rate_perm1 = len(IDs210min_inBaseline)
perc_multifibre_bundles_210min_rate_perm1 = nr_multifibre_bundles_210min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_210_perm2 = np.arange(30,1291,210)
data_multi_210min_rate = data_multi[data_multi['Time'].isin(times_210_perm2)] 
unique_multifibres_210min_rate  = data_multi_210min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_210min = unique_multifibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_multifibres_210min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_210min_rate_perm2 = len(IDs210min_inBaseline)
perc_multifibre_bundles_210min_rate_perm2 = nr_multifibre_bundles_210min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_210_perm3 = np.arange(60,1291,210)
data_multi_210min_rate = data_multi[data_multi['Time'].isin(times_210_perm3)] 
unique_multifibres_210min_rate  = data_multi_210min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_210min = unique_multifibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_multifibres_210min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_210min_rate_perm3 = len(IDs210min_inBaseline)
perc_multifibre_bundles_210min_rate_perm3 = nr_multifibre_bundles_210min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_210_perm4 = np.arange(90,1291,210)
data_multi_210min_rate = data_multi[data_multi['Time'].isin(times_210_perm4)] 
unique_multifibres_210min_rate  = data_multi_210min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_210min = unique_multifibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_multifibres_210min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_210min_rate_perm4 = len(IDs210min_inBaseline)
perc_multifibre_bundles_210min_rate_perm4 = nr_multifibre_bundles_210min_rate_perm4 / total_nr_multifibre_bundles

#permuation 5
times_210_perm5 = np.arange(120,1291,210)
data_multi_210min_rate = data_multi[data_multi['Time'].isin(times_210_perm5)] 
unique_multifibres_210min_rate  = data_multi_210min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_210min = unique_multifibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_multifibres_210min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_210min_rate_perm5 = len(IDs210min_inBaseline)
perc_multifibre_bundles_210min_rate_perm5 = nr_multifibre_bundles_210min_rate_perm5 / total_nr_multifibre_bundles

#permuation 6
times_210_perm6 = np.arange(150,1291,210)
data_multi_210min_rate = data_multi[data_multi['Time'].isin(times_210_perm6)] 
unique_multifibres_210min_rate  = data_multi_210min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_210min = unique_multifibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_multifibres_210min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_210min_rate_perm6 = len(IDs210min_inBaseline)
perc_multifibre_bundles_210min_rate_perm6 = nr_multifibre_bundles_210min_rate_perm6 / total_nr_multifibre_bundles

#permuation 7
times_210_perm7 = np.arange(180,1291,210)
data_multi_210min_rate = data_multi[data_multi['Time'].isin(times_210_perm7)] 
unique_multifibres_210min_rate  = data_multi_210min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_210min = unique_multifibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_multifibres_210min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_210min_rate_perm7 = len(IDs210min_inBaseline)
perc_multifibre_bundles_210min_rate_perm7 = nr_multifibre_bundles_210min_rate_perm7 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_210min_rate = np.mean([nr_multifibre_bundles_210min_rate_perm1, nr_multifibre_bundles_210min_rate_perm2, 
                                            nr_multifibre_bundles_210min_rate_perm3, nr_multifibre_bundles_210min_rate_perm4,
                                            nr_multifibre_bundles_210min_rate_perm5, nr_multifibre_bundles_210min_rate_perm6,
                                            nr_multifibre_bundles_210min_rate_perm7])
perc_multifibre_bundles_210min_rate = np.mean([perc_multifibre_bundles_210min_rate_perm1, perc_multifibre_bundles_210min_rate_perm2, 
                                              perc_multifibre_bundles_210min_rate_perm3, perc_multifibre_bundles_210min_rate_perm4,
                                              perc_multifibre_bundles_210min_rate_perm5, perc_multifibre_bundles_210min_rate_perm6,
                                              perc_multifibre_bundles_210min_rate_perm7])


frame_hour_rates = [30,60,90,120,150, 180, 210]
nr_detected_multifibre_bundles_perm = [total_nr_multifibre_bundles,
                                  nr_multifibre_bundles_60min_rate, nr_multifibre_bundles_90min_rate,
                                  nr_multifibre_bundles_120min_rate, nr_multifibre_bundles_150min_rate,
                                  nr_multifibre_bundles_180min_rate, nr_multifibre_bundles_210min_rate]
percentage_detected_multifibre_bundles_perm = [perc_total_nr_multifibre_bundles, perc_multifibre_bundles_60min_rate,
                                          perc_multifibre_bundles_90min_rate, perc_multifibre_bundles_120min_rate,
                                          perc_multifibre_bundles_150min_rate, perc_multifibre_bundles_180min_rate,
                                          perc_multifibre_bundles_210min_rate]


data_frame_rates_perm = pd.DataFrame(list(zip(frame_hour_rates, nr_detected_multifibre_bundles_perm, percentage_detected_multifibre_bundles_perm)),
                                columns = ['Frame Rate', 'Nr detected multifibre bundles', 'Percentage detected multifibre bundles'])

x = data_frame_rates_perm['Frame Rate']
y = data_frame_rates_perm['Percentage detected multifibre bundles']
popt, pcov = curve_fit(func, x, y,  p0 = [1, 80, 0.2, 1], maxfev=50000)

# # calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 0 min
# renorm_factor = func(0, *popt)

# calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 5 min
renorm_factor = func(5, *popt)

x_continuous = np.arange(30, 220, 1)
# Plot the frame rate and the number of detected multifibre bundles

x_continuous = np.arange(0, 220, 1)

plt.figure(figsize=(9, 7))
fs =18

width = 0.3

plt.rc('grid', linestyle="-", color='lightgrey')
plt.grid(True)
plt.scatter(data_frame_rates_perm['Frame Rate'], data_frame_rates_perm['Percentage detected multifibre bundles'], color = 'black', label='Data')
plt.plot(x_continuous, func(x_continuous, *popt), label='Inverse logistic regression', color = 'tomato')
plt.legend(fontsize=14)


plt.xlabel('Frame Rate (minutes)', fontsize=fs)
plt.ylabel('Proportion of detected multi-fibre axonal structures', fontsize=fs)
#plt.title('Detected multifibre axonal structures for different frame rates', fontsize = 18)

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.set_xlim([0, 220])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,30,60,90,120,150,180,210])
plt.show(block=True)
 
plt.close()

# Adjusted to renormalization factor

x_continuous = np.arange(0, 220, 1)

plt.figure(figsize=(9, 7))


fs =18

width = 0.3

plt.rc('grid', linestyle="-", color='lightgrey')
plt.grid(True)
plt.scatter(data_frame_rates_perm['Frame Rate'], data_frame_rates_perm['Percentage detected multifibre bundles']/renorm_factor, color = 'black', label='Data, renormalized')
plt.plot(x_continuous, func(x_continuous, *popt)/renorm_factor, label='Inverse logistic regression', color = 'tomato')
plt.legend(fontsize=14)


plt.xlabel('Frame Rate (minutes)', fontsize=fs)
plt.ylabel('Proportion of detected multi-fibre axonal structures', fontsize=fs)
#plt.title('Detected multifibre axonal structures for different frame rates', fontsize = 18)

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.set_xlim([0, 220])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,30,60,90,120,150,180,210])

#ax.set_xticklabels([0,5,20,40,60,80,100,120,140])
plt.show(block=True)
 
plt.close()

perc_detected_30min = func(30, *popt)/renorm_factor * 100
perc_detected_60min = func(60, *popt)/renorm_factor * 100
difference_perc_detected_30_60 = perc_detected_30min - perc_detected_60min 

print('Percentage of multi-fibre axonal structures being detected with a 60 minute frame rate: ' + str(perc_detected_60min) + '%')
print('Percentage of multi-fibre axonal structures being detected with a 30 minute frame rate: ' + str(perc_detected_30min) + '%')
print('Difference of percentage of multi-fibre axonal structures being detected between 30 and 60 minute frame rate: ' + str(difference_perc_detected_30_60) + '%')



# %% singlefibre axonal structures

# Times are now from 0 to 1290 minutes, with a 30 minute frame rate


########### Considering permutations #######################

# data we get if we measure every 60 minutes

#permutation 1
times_60_perm1 = np.arange(0,1291,60)
data_single_60min_rate = data_single[data_single['Time'].isin(times_60_perm1)] 
unique_singlefibres_60min_rate  = data_single_60min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_60min = unique_singlefibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_singlefibres_60min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_60min_rate_perm1 = len(IDs60min_inBaseline)
perc_singlefibre_bundles_60min_rate_perm1 = nr_singlefibre_bundles_60min_rate_perm1 / total_nr_singlefibre_bundles

#permutation 2
times_60_perm2 = np.arange(30,1291,60)
data_single_60min_rate = data_single[data_single['Time'].isin(times_60_perm2)] 
unique_singlefibres_60min_rate  = data_single_60min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_60min = unique_singlefibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_singlefibres_60min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_60min_rate_perm2 = len(IDs60min_inBaseline)
perc_singlefibre_bundles_60min_rate_perm2 = nr_singlefibre_bundles_60min_rate_perm2 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_60min_rate = np.mean([nr_singlefibre_bundles_60min_rate_perm1, nr_singlefibre_bundles_60min_rate_perm2])
perc_singlefibre_bundles_60min_rate = np.mean([perc_singlefibre_bundles_60min_rate_perm1, perc_singlefibre_bundles_60min_rate_perm2])



# data we get if we measure every 90 minutes

#permutation 1
times_90_perm1 = np.arange(0,1291,90)
data_single_90min_rate = data_single[data_single['Time'].isin(times_90_perm1)] 
unique_singlefibres_90min_rate  = data_single_90min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_90min = unique_singlefibres_90min_rate['TrackID']
IDs90min_inBaseline = list_trackID_singlefibres_90min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_90min_rate_perm1 = len(IDs90min_inBaseline)
perc_singlefibre_bundles_90min_rate_perm1 = nr_singlefibre_bundles_90min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_90_perm2 = np.arange(30,1291,90)
data_single_90min_rate = data_single[data_single['Time'].isin(times_90_perm2)] 
unique_singlefibres_90min_rate  = data_single_90min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_90min = unique_singlefibres_90min_rate['TrackID']
IDs90min_inBaseline = list_trackID_singlefibres_90min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_90min_rate_perm2 = len(IDs90min_inBaseline)
perc_singlefibre_bundles_90min_rate_perm2 = nr_singlefibre_bundles_90min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_90_perm3 = np.arange(60,1291,90)
data_single_90min_rate = data_single[data_single['Time'].isin(times_90_perm3)] 
unique_singlefibres_90min_rate  = data_single_90min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_90min = unique_singlefibres_90min_rate['TrackID']
IDs90min_inBaseline = list_trackID_singlefibres_90min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_90min_rate_perm3 = len(IDs90min_inBaseline)
perc_singlefibre_bundles_90min_rate_perm3 = nr_singlefibre_bundles_90min_rate_perm3 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_90min_rate = np.mean([nr_singlefibre_bundles_90min_rate_perm1, nr_singlefibre_bundles_90min_rate_perm2, 
                                            nr_singlefibre_bundles_90min_rate_perm3])
perc_singlefibre_bundles_90min_rate = np.mean([perc_singlefibre_bundles_90min_rate_perm1, perc_singlefibre_bundles_90min_rate_perm2, 
                                              perc_singlefibre_bundles_90min_rate_perm3])



# data we get if we measure every 120 minutes

#permutation 1
times_120_perm1 = np.arange(0,1291,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm1)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm1 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm1 = nr_singlefibre_bundles_120min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_120_perm2 = np.arange(30,1291,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm2)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm2 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm2 = nr_singlefibre_bundles_120min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_120_perm3 = np.arange(60,1291,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm3)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm3 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm3 = nr_singlefibre_bundles_120min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_120_perm4 = np.arange(90,1291,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm4)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm4 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm4 = nr_singlefibre_bundles_120min_rate_perm4 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_120min_rate = np.mean([nr_singlefibre_bundles_120min_rate_perm1, nr_singlefibre_bundles_120min_rate_perm2, 
                                            nr_singlefibre_bundles_120min_rate_perm3, nr_singlefibre_bundles_120min_rate_perm4])
perc_singlefibre_bundles_120min_rate = np.mean([perc_singlefibre_bundles_120min_rate_perm1, perc_singlefibre_bundles_120min_rate_perm2, 
                                              perc_singlefibre_bundles_120min_rate_perm3, perc_singlefibre_bundles_120min_rate_perm4])


# data we get if we measure every 150 minutes

#permutation 1
times_150_perm1 = np.arange(0,1291,150)
data_single_150min_rate = data_single[data_single['Time'].isin(times_150_perm1)] 
unique_singlefibres_150min_rate  = data_single_150min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_150min = unique_singlefibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_singlefibres_150min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_150min_rate_perm1 = len(IDs150min_inBaseline)
perc_singlefibre_bundles_150min_rate_perm1 = nr_singlefibre_bundles_150min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_150_perm2 = np.arange(30,1291,150)
data_single_150min_rate = data_single[data_single['Time'].isin(times_150_perm2)] 
unique_singlefibres_150min_rate  = data_single_150min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_150min = unique_singlefibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_singlefibres_150min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_150min_rate_perm2 = len(IDs150min_inBaseline)
perc_singlefibre_bundles_150min_rate_perm2 = nr_singlefibre_bundles_150min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_150_perm3 = np.arange(60,1291,150)
data_single_150min_rate = data_single[data_single['Time'].isin(times_150_perm3)] 
unique_singlefibres_150min_rate  = data_single_150min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_150min = unique_singlefibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_singlefibres_150min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_150min_rate_perm3 = len(IDs150min_inBaseline)
perc_singlefibre_bundles_150min_rate_perm3 = nr_singlefibre_bundles_150min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_150_perm4 = np.arange(90,1291,150)
data_single_150min_rate = data_single[data_single['Time'].isin(times_150_perm4)] 
unique_singlefibres_150min_rate  = data_single_150min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_150min = unique_singlefibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_singlefibres_150min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_150min_rate_perm4 = len(IDs150min_inBaseline)
perc_singlefibre_bundles_150min_rate_perm4 = nr_singlefibre_bundles_150min_rate_perm4 / total_nr_singlefibre_bundles

#permuation 5
times_150_perm5 = np.arange(120,1291,150)
data_single_150min_rate = data_single[data_single['Time'].isin(times_150_perm5)] 
unique_singlefibres_150min_rate  = data_single_150min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_150min = unique_singlefibres_150min_rate['TrackID']
IDs150min_inBaseline = list_trackID_singlefibres_150min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_150min_rate_perm5 = len(IDs150min_inBaseline)
perc_singlefibre_bundles_150min_rate_perm5 = nr_singlefibre_bundles_150min_rate_perm5 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_150min_rate = np.mean([nr_singlefibre_bundles_150min_rate_perm1, nr_singlefibre_bundles_150min_rate_perm2, 
                                            nr_singlefibre_bundles_150min_rate_perm3, nr_singlefibre_bundles_150min_rate_perm4,
                                            nr_singlefibre_bundles_150min_rate_perm5])
perc_singlefibre_bundles_150min_rate = np.mean([perc_singlefibre_bundles_150min_rate_perm1, perc_singlefibre_bundles_150min_rate_perm2, 
                                              perc_singlefibre_bundles_150min_rate_perm3, perc_singlefibre_bundles_150min_rate_perm4,
                                              perc_singlefibre_bundles_150min_rate_perm5])

# data we get if we measure every 180 minutes

#permutation 1
times_180_perm1 = np.arange(0,1291,180)
data_single_180min_rate = data_single[data_single['Time'].isin(times_180_perm1)] 
unique_singlefibres_180min_rate  = data_single_180min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_180min = unique_singlefibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_singlefibres_180min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_180min_rate_perm1 = len(IDs180min_inBaseline)
perc_singlefibre_bundles_180min_rate_perm1 = nr_singlefibre_bundles_180min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_180_perm2 = np.arange(30,1291,180)
data_single_180min_rate = data_single[data_single['Time'].isin(times_180_perm2)] 
unique_singlefibres_180min_rate  = data_single_180min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_180min = unique_singlefibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_singlefibres_180min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_180min_rate_perm2 = len(IDs180min_inBaseline)
perc_singlefibre_bundles_180min_rate_perm2 = nr_singlefibre_bundles_180min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_180_perm3 = np.arange(60,1291,180)
data_single_180min_rate = data_single[data_single['Time'].isin(times_180_perm3)] 
unique_singlefibres_180min_rate  = data_single_180min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_180min = unique_singlefibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_singlefibres_180min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_180min_rate_perm3 = len(IDs180min_inBaseline)
perc_singlefibre_bundles_180min_rate_perm3 = nr_singlefibre_bundles_180min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_180_perm4 = np.arange(90,1291,180)
data_single_180min_rate = data_single[data_single['Time'].isin(times_180_perm4)] 
unique_singlefibres_180min_rate  = data_single_180min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_180min = unique_singlefibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_singlefibres_180min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_180min_rate_perm4 = len(IDs180min_inBaseline)
perc_singlefibre_bundles_180min_rate_perm4 = nr_singlefibre_bundles_180min_rate_perm4 / total_nr_singlefibre_bundles

#permuation 5
times_180_perm5 = np.arange(120,1291,180)
data_single_180min_rate = data_single[data_single['Time'].isin(times_180_perm5)] 
unique_singlefibres_180min_rate  = data_single_180min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_180min = unique_singlefibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_singlefibres_180min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_180min_rate_perm5 = len(IDs180min_inBaseline)
perc_singlefibre_bundles_180min_rate_perm5 = nr_singlefibre_bundles_180min_rate_perm5 / total_nr_singlefibre_bundles

#permuation 6
times_180_perm6 = np.arange(150,1291,180)
data_single_180min_rate = data_single[data_single['Time'].isin(times_180_perm6)] 
unique_singlefibres_180min_rate  = data_single_180min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_180min = unique_singlefibres_180min_rate['TrackID']
IDs180min_inBaseline = list_trackID_singlefibres_180min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_180min_rate_perm6 = len(IDs180min_inBaseline)
perc_singlefibre_bundles_180min_rate_perm6 = nr_singlefibre_bundles_180min_rate_perm6 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_180min_rate = np.mean([nr_singlefibre_bundles_180min_rate_perm1, nr_singlefibre_bundles_180min_rate_perm2, 
                                            nr_singlefibre_bundles_180min_rate_perm3, nr_singlefibre_bundles_180min_rate_perm4,
                                            nr_singlefibre_bundles_180min_rate_perm5, nr_singlefibre_bundles_180min_rate_perm6])
perc_singlefibre_bundles_180min_rate = np.mean([perc_singlefibre_bundles_180min_rate_perm1, perc_singlefibre_bundles_180min_rate_perm2, 
                                              perc_singlefibre_bundles_180min_rate_perm3, perc_singlefibre_bundles_180min_rate_perm4,
                                              perc_singlefibre_bundles_180min_rate_perm5, perc_singlefibre_bundles_180min_rate_perm6])


# data we get if we measure every 210 minutes

#permutation 1
times_210_perm1 = np.arange(0,1291,210)
data_single_210min_rate = data_single[data_single['Time'].isin(times_210_perm1)] 
unique_singlefibres_210min_rate  = data_single_210min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_210min = unique_singlefibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_singlefibres_210min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_210min_rate_perm1 = len(IDs210min_inBaseline)
perc_singlefibre_bundles_210min_rate_perm1 = nr_singlefibre_bundles_210min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_210_perm2 = np.arange(30,1291,210)
data_single_210min_rate = data_single[data_single['Time'].isin(times_210_perm2)] 
unique_singlefibres_210min_rate  = data_single_210min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_210min = unique_singlefibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_singlefibres_210min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_210min_rate_perm2 = len(IDs210min_inBaseline)
perc_singlefibre_bundles_210min_rate_perm2 = nr_singlefibre_bundles_210min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_210_perm3 = np.arange(60,1291,210)
data_single_210min_rate = data_single[data_single['Time'].isin(times_210_perm3)] 
unique_singlefibres_210min_rate  = data_single_210min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_210min = unique_singlefibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_singlefibres_210min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_210min_rate_perm3 = len(IDs210min_inBaseline)
perc_singlefibre_bundles_210min_rate_perm3 = nr_singlefibre_bundles_210min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_210_perm4 = np.arange(90,1291,210)
data_single_210min_rate = data_single[data_single['Time'].isin(times_210_perm4)] 
unique_singlefibres_210min_rate  = data_single_210min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_210min = unique_singlefibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_singlefibres_210min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_210min_rate_perm4 = len(IDs210min_inBaseline)
perc_singlefibre_bundles_210min_rate_perm4 = nr_singlefibre_bundles_210min_rate_perm4 / total_nr_singlefibre_bundles

#permuation 5
times_210_perm5 = np.arange(120,1291,210)
data_single_210min_rate = data_single[data_single['Time'].isin(times_210_perm5)] 
unique_singlefibres_210min_rate  = data_single_210min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_210min = unique_singlefibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_singlefibres_210min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_210min_rate_perm5 = len(IDs210min_inBaseline)
perc_singlefibre_bundles_210min_rate_perm5 = nr_singlefibre_bundles_210min_rate_perm5 / total_nr_singlefibre_bundles

#permuation 6
times_210_perm6 = np.arange(150,1291,210)
data_single_210min_rate = data_single[data_single['Time'].isin(times_210_perm6)] 
unique_singlefibres_210min_rate  = data_single_210min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_210min = unique_singlefibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_singlefibres_210min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_210min_rate_perm6 = len(IDs210min_inBaseline)
perc_singlefibre_bundles_210min_rate_perm6 = nr_singlefibre_bundles_210min_rate_perm6 / total_nr_singlefibre_bundles

#permuation 7
times_210_perm7 = np.arange(180,1291,210)
data_single_210min_rate = data_single[data_single['Time'].isin(times_210_perm7)] 
unique_singlefibres_210min_rate  = data_single_210min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_210min = unique_singlefibres_210min_rate['TrackID']
IDs210min_inBaseline = list_trackID_singlefibres_210min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_210min_rate_perm7 = len(IDs210min_inBaseline)
perc_singlefibre_bundles_210min_rate_perm7 = nr_singlefibre_bundles_210min_rate_perm7 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_210min_rate = np.mean([nr_singlefibre_bundles_210min_rate_perm1, nr_singlefibre_bundles_210min_rate_perm2, 
                                            nr_singlefibre_bundles_210min_rate_perm3, nr_singlefibre_bundles_210min_rate_perm4,
                                            nr_singlefibre_bundles_210min_rate_perm5, nr_singlefibre_bundles_210min_rate_perm6,
                                            nr_singlefibre_bundles_210min_rate_perm7])
perc_singlefibre_bundles_210min_rate = np.mean([perc_singlefibre_bundles_210min_rate_perm1, perc_singlefibre_bundles_210min_rate_perm2, 
                                              perc_singlefibre_bundles_210min_rate_perm3, perc_singlefibre_bundles_210min_rate_perm4,
                                              perc_singlefibre_bundles_210min_rate_perm5, perc_singlefibre_bundles_210min_rate_perm6,
                                              perc_singlefibre_bundles_210min_rate_perm7])


frame_hour_rates = [30,60,90,120,150, 180, 210]
nr_detected_singlefibre_bundles_perm = [total_nr_singlefibre_bundles,
                                  nr_singlefibre_bundles_60min_rate, nr_singlefibre_bundles_90min_rate,
                                  nr_singlefibre_bundles_120min_rate, nr_singlefibre_bundles_150min_rate,
                                  nr_singlefibre_bundles_180min_rate, nr_singlefibre_bundles_210min_rate]
percentage_detected_singlefibre_bundles_perm = [perc_total_nr_singlefibre_bundles, perc_singlefibre_bundles_60min_rate,
                                          perc_singlefibre_bundles_90min_rate, perc_singlefibre_bundles_120min_rate,
                                          perc_singlefibre_bundles_150min_rate, perc_singlefibre_bundles_180min_rate,
                                          perc_singlefibre_bundles_210min_rate]


data_frame_rates_perm = pd.DataFrame(list(zip(frame_hour_rates, nr_detected_singlefibre_bundles_perm, percentage_detected_singlefibre_bundles_perm)),
                                columns = ['Frame Rate', 'Nr detected singlefibre bundles', 'Percentage detected singlefibre bundles'])

x = data_frame_rates_perm['Frame Rate']
y = data_frame_rates_perm['Percentage detected singlefibre bundles']
popt, pcov = curve_fit(func, x, y,  p0 = [1, 80, 0.2, 1], maxfev=50000)

# # calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 0 min
# renorm_factor = func(0, *popt)

# calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 5 min
renorm_factor = func(5, *popt)

x_continuous = np.arange(30, 220, 1)
# Plot the frame rate and the number of detected singlefibre bundles

x_continuous = np.arange(0, 220, 1)

plt.figure(figsize=(9, 7))
fs =18

width = 0.3

plt.rc('grid', linestyle="-", color='lightgrey')
plt.grid(True)
plt.scatter(data_frame_rates_perm['Frame Rate'], data_frame_rates_perm['Percentage detected singlefibre bundles'], color = 'black', label='Data')
plt.plot(x_continuous, func(x_continuous, *popt), label='Inverse logistic regression', color = 'tomato')
plt.legend(fontsize=14)


plt.xlabel('Frame Rate (minutes)', fontsize=fs)
plt.ylabel('Proportion of detected single-fibre axonal structures', fontsize=fs)
#plt.title('Detected singlefibre axonal structures for different frame rates', fontsize = 18)

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.set_xlim([0, 220])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,30,60,90,120,150,180,210])

plt.show(block=True)
 
plt.close()

# Adjusted to renormalization factor

x_continuous = np.arange(0, 220, 1)

plt.figure(figsize=(9, 7))


fs =18

width = 0.3

plt.rc('grid', linestyle="-", color='lightgrey')
plt.grid(True)
plt.scatter(data_frame_rates_perm['Frame Rate'], data_frame_rates_perm['Percentage detected singlefibre bundles']/renorm_factor, color = 'black', label='Data, renormalized')
plt.plot(x_continuous, func(x_continuous, *popt)/renorm_factor, label='Inverse logistic regression', color = 'tomato')
plt.legend(fontsize=14)


plt.xlabel('Frame Rate (minutes)', fontsize=fs)
plt.ylabel('Proportion of detected single-fibre axonal structures', fontsize=fs)
#plt.title('Detected singlefibre axonal structures for different frame rates', fontsize = 18)

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.set_xlim([0, 220])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,30,60,90,120,150,180,210])

#ax.set_xticklabels([0,5,20,40,60,80,100,120,140])
plt.show(block=True)
 
plt.close()


