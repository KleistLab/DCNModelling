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
filename = os.path.join(dirname, 'dataset-20231211.csv')

data = pd.read_csv(filename, sep = ';')
data = data.iloc[:, :-1]


# Remove reference lines from the data
data = data[data.Name == 'structure']

# TEST- keep only multifibre structures- option 1
data_multi = data[data.Position_Type_Track == 'multi']

data_single = data[data.Position_Type_Track == 'single']


# TEST- keep only multifibre structures- option 2
# data['Multifibre'] = data['Fibers'].apply(lambda x: x>1)
# data['Multifibre'] = data['TrackID'].map(data.groupby('TrackID')['Multifibre'].any())

# data_multi2 = data[data.Multifibre == True]

# CONCLUSION: both options to get multifibres give same results, so we keep the
# simplest one

# 
data_multi_unique  = data_multi.drop_duplicates('TrackID')

data_single_unique  = data_single.drop_duplicates('TrackID')


# Total number of axonal structuress that had more than one fibre (multifibre) at some point in time
total_nr_multifibre_bundles = len(data_multi_unique)
list_trackID_multifibres_baseline = data_multi_unique['TrackID']

perc_total_nr_multifibre_bundles = total_nr_multifibre_bundles / total_nr_multifibre_bundles

# Total number of axonal structuress that had more only one fibre (single-fibre) at all time points
total_nr_singlefibre_bundles = len(data_single_unique)
list_trackID_singlefibres_baseline = data_single_unique['TrackID']

perc_total_nr_singlefibre_bundles = total_nr_singlefibre_bundles / total_nr_singlefibre_bundles


# %% Multifibre axonal structures

# Times are now from 0 to 600 minutes, with a 20 minute frame rate


########### Considering permutations #######################

# data we get if we measure every 40 minutes

#permutation 1
times_40_perm1 = np.arange(0,601,40)
data_multi_40min_rate = data_multi[data_multi['Time'].isin(times_40_perm1)] 
unique_multifibres_40min_rate  = data_multi_40min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_40min = unique_multifibres_40min_rate['TrackID']
IDs40min_inBaseline = list_trackID_multifibres_40min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_40min_rate_perm1 = len(IDs40min_inBaseline)
perc_multifibre_bundles_40min_rate_perm1 = nr_multifibre_bundles_40min_rate_perm1 / total_nr_multifibre_bundles

#permutation 2
times_40_perm2 = np.arange(20,601,40)
data_multi_40min_rate = data_multi[data_multi['Time'].isin(times_40_perm2)] 
unique_multifibres_40min_rate  = data_multi_40min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_40min = unique_multifibres_40min_rate['TrackID']
IDs40min_inBaseline = list_trackID_multifibres_40min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_40min_rate_perm2 = len(IDs40min_inBaseline)
perc_multifibre_bundles_40min_rate_perm2 = nr_multifibre_bundles_40min_rate_perm2 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_40min_rate = np.mean([nr_multifibre_bundles_40min_rate_perm1, nr_multifibre_bundles_40min_rate_perm2])
perc_multifibre_bundles_40min_rate = np.mean([perc_multifibre_bundles_40min_rate_perm1, perc_multifibre_bundles_40min_rate_perm2])


# data we get if we measure every 60 minutes

#permutation 1
times_60_perm1 = np.arange(0,601,60)
data_multi_60min_rate = data_multi[data_multi['Time'].isin(times_60_perm1)] 
unique_multifibres_60min_rate  = data_multi_60min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_60min = unique_multifibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_multifibres_60min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_60min_rate_perm1 = len(IDs60min_inBaseline)
perc_multifibre_bundles_60min_rate_perm1 = nr_multifibre_bundles_60min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_60_perm2 = np.arange(20,601,60)
data_multi_60min_rate = data_multi[data_multi['Time'].isin(times_60_perm2)] 
unique_multifibres_60min_rate  = data_multi_60min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_60min = unique_multifibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_multifibres_60min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_60min_rate_perm2 = len(IDs60min_inBaseline)
perc_multifibre_bundles_60min_rate_perm2 = nr_multifibre_bundles_60min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_60_perm3 = np.arange(40,601,60)
data_multi_60min_rate = data_multi[data_multi['Time'].isin(times_60_perm3)] 
unique_multifibres_60min_rate  = data_multi_60min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_60min = unique_multifibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_multifibres_60min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_60min_rate_perm3 = len(IDs60min_inBaseline)
perc_multifibre_bundles_60min_rate_perm3 = nr_multifibre_bundles_60min_rate_perm3 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_60min_rate = np.mean([nr_multifibre_bundles_60min_rate_perm1, nr_multifibre_bundles_60min_rate_perm2, 
                                            nr_multifibre_bundles_60min_rate_perm3])
perc_multifibre_bundles_60min_rate = np.mean([perc_multifibre_bundles_60min_rate_perm1, perc_multifibre_bundles_60min_rate_perm2, 
                                              perc_multifibre_bundles_60min_rate_perm3])



# data we get if we measure every 80 minutes

#permutation 1
times_80_perm1 = np.arange(0,601,80)
data_multi_80min_rate = data_multi[data_multi['Time'].isin(times_80_perm1)] 
unique_multifibres_80min_rate  = data_multi_80min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_80min = unique_multifibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_multifibres_80min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_80min_rate_perm1 = len(IDs80min_inBaseline)
perc_multifibre_bundles_80min_rate_perm1 = nr_multifibre_bundles_80min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_80_perm2 = np.arange(20,601,80)
data_multi_80min_rate = data_multi[data_multi['Time'].isin(times_80_perm2)] 
unique_multifibres_80min_rate  = data_multi_80min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_80min = unique_multifibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_multifibres_80min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_80min_rate_perm2 = len(IDs80min_inBaseline)
perc_multifibre_bundles_80min_rate_perm2 = nr_multifibre_bundles_80min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_80_perm3 = np.arange(40,601,80)
data_multi_80min_rate = data_multi[data_multi['Time'].isin(times_80_perm3)] 
unique_multifibres_80min_rate  = data_multi_80min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_80min = unique_multifibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_multifibres_80min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_80min_rate_perm3 = len(IDs80min_inBaseline)
perc_multifibre_bundles_80min_rate_perm3 = nr_multifibre_bundles_80min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_80_perm4 = np.arange(60,601,80)
data_multi_80min_rate = data_multi[data_multi['Time'].isin(times_80_perm4)] 
unique_multifibres_80min_rate  = data_multi_80min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_80min = unique_multifibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_multifibres_80min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_80min_rate_perm4 = len(IDs80min_inBaseline)
perc_multifibre_bundles_80min_rate_perm4 = nr_multifibre_bundles_80min_rate_perm4 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_80min_rate = np.mean([nr_multifibre_bundles_80min_rate_perm1, nr_multifibre_bundles_80min_rate_perm2, 
                                            nr_multifibre_bundles_80min_rate_perm3, nr_multifibre_bundles_80min_rate_perm4])
perc_multifibre_bundles_80min_rate = np.mean([perc_multifibre_bundles_80min_rate_perm1, perc_multifibre_bundles_80min_rate_perm2, 
                                              perc_multifibre_bundles_80min_rate_perm3, perc_multifibre_bundles_80min_rate_perm4])


# data we get if we measure every 100 minutes

#permutation 1
times_100_perm1 = np.arange(0,601,100)
data_multi_100min_rate = data_multi[data_multi['Time'].isin(times_100_perm1)] 
unique_multifibres_100min_rate  = data_multi_100min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_100min = unique_multifibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_multifibres_100min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_100min_rate_perm1 = len(IDs100min_inBaseline)
perc_multifibre_bundles_100min_rate_perm1 = nr_multifibre_bundles_100min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_100_perm2 = np.arange(20,601,100)
data_multi_100min_rate = data_multi[data_multi['Time'].isin(times_100_perm2)] 
unique_multifibres_100min_rate  = data_multi_100min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_100min = unique_multifibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_multifibres_100min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_100min_rate_perm2 = len(IDs100min_inBaseline)
perc_multifibre_bundles_100min_rate_perm2 = nr_multifibre_bundles_100min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_100_perm3 = np.arange(40,601,100)
data_multi_100min_rate = data_multi[data_multi['Time'].isin(times_100_perm3)] 
unique_multifibres_100min_rate  = data_multi_100min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_100min = unique_multifibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_multifibres_100min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_100min_rate_perm3 = len(IDs100min_inBaseline)
perc_multifibre_bundles_100min_rate_perm3 = nr_multifibre_bundles_100min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_100_perm4 = np.arange(60,601,100)
data_multi_100min_rate = data_multi[data_multi['Time'].isin(times_100_perm4)] 
unique_multifibres_100min_rate  = data_multi_100min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_100min = unique_multifibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_multifibres_100min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_100min_rate_perm4 = len(IDs100min_inBaseline)
perc_multifibre_bundles_100min_rate_perm4 = nr_multifibre_bundles_100min_rate_perm4 / total_nr_multifibre_bundles

#permuation 5
times_100_perm5 = np.arange(80,601,100)
data_multi_100min_rate = data_multi[data_multi['Time'].isin(times_100_perm5)] 
unique_multifibres_100min_rate  = data_multi_100min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_100min = unique_multifibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_multifibres_100min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_100min_rate_perm5 = len(IDs100min_inBaseline)
perc_multifibre_bundles_100min_rate_perm5 = nr_multifibre_bundles_100min_rate_perm5 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_100min_rate = np.mean([nr_multifibre_bundles_100min_rate_perm1, nr_multifibre_bundles_100min_rate_perm2, 
                                            nr_multifibre_bundles_100min_rate_perm3, nr_multifibre_bundles_100min_rate_perm4,
                                            nr_multifibre_bundles_100min_rate_perm5])
perc_multifibre_bundles_100min_rate = np.mean([perc_multifibre_bundles_100min_rate_perm1, perc_multifibre_bundles_100min_rate_perm2, 
                                              perc_multifibre_bundles_100min_rate_perm3, perc_multifibre_bundles_100min_rate_perm4,
                                              perc_multifibre_bundles_100min_rate_perm5])

# data we get if we measure every 120 minutes

#permutation 1
times_120_perm1 = np.arange(0,601,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm1)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm1 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm1 = nr_multifibre_bundles_120min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_120_perm2 = np.arange(20,601,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm2)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm2 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm2 = nr_multifibre_bundles_120min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_120_perm3 = np.arange(40,601,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm3)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm3 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm3 = nr_multifibre_bundles_120min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_120_perm4 = np.arange(60,601,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm4)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm4 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm4 = nr_multifibre_bundles_120min_rate_perm4 / total_nr_multifibre_bundles

#permuation 5
times_120_perm5 = np.arange(80,601,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm5)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm5 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm5 = nr_multifibre_bundles_120min_rate_perm5 / total_nr_multifibre_bundles

#permuation 6
times_120_perm6 = np.arange(100,601,120)
data_multi_120min_rate = data_multi[data_multi['Time'].isin(times_120_perm6)] 
unique_multifibres_120min_rate  = data_multi_120min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_120min = unique_multifibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_multifibres_120min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_120min_rate_perm6 = len(IDs120min_inBaseline)
perc_multifibre_bundles_120min_rate_perm6 = nr_multifibre_bundles_120min_rate_perm6 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_120min_rate = np.mean([nr_multifibre_bundles_120min_rate_perm1, nr_multifibre_bundles_120min_rate_perm2, 
                                            nr_multifibre_bundles_120min_rate_perm3, nr_multifibre_bundles_120min_rate_perm4,
                                            nr_multifibre_bundles_120min_rate_perm5, nr_multifibre_bundles_120min_rate_perm6])
perc_multifibre_bundles_120min_rate = np.mean([perc_multifibre_bundles_120min_rate_perm1, perc_multifibre_bundles_120min_rate_perm2, 
                                              perc_multifibre_bundles_120min_rate_perm3, perc_multifibre_bundles_120min_rate_perm4,
                                              perc_multifibre_bundles_120min_rate_perm5, perc_multifibre_bundles_120min_rate_perm6])


# data we get if we measure every 140 minutes

#permutation 1
times_140_perm1 = np.arange(0,601,140)
data_multi_140min_rate = data_multi[data_multi['Time'].isin(times_140_perm1)] 
unique_multifibres_140min_rate  = data_multi_140min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_140min = unique_multifibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_multifibres_140min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_140min_rate_perm1 = len(IDs140min_inBaseline)
perc_multifibre_bundles_140min_rate_perm1 = nr_multifibre_bundles_140min_rate_perm1 / total_nr_multifibre_bundles

#permuation 2
times_140_perm2 = np.arange(20,601,140)
data_multi_140min_rate = data_multi[data_multi['Time'].isin(times_140_perm2)] 
unique_multifibres_140min_rate  = data_multi_140min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_140min = unique_multifibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_multifibres_140min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_140min_rate_perm2 = len(IDs140min_inBaseline)
perc_multifibre_bundles_140min_rate_perm2 = nr_multifibre_bundles_140min_rate_perm2 / total_nr_multifibre_bundles

#permuation 3
times_140_perm3 = np.arange(40,601,140)
data_multi_140min_rate = data_multi[data_multi['Time'].isin(times_140_perm3)] 
unique_multifibres_140min_rate  = data_multi_140min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_140min = unique_multifibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_multifibres_140min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_140min_rate_perm3 = len(IDs140min_inBaseline)
perc_multifibre_bundles_140min_rate_perm3 = nr_multifibre_bundles_140min_rate_perm3 / total_nr_multifibre_bundles

#permuation 4
times_140_perm4 = np.arange(60,601,140)
data_multi_140min_rate = data_multi[data_multi['Time'].isin(times_140_perm4)] 
unique_multifibres_140min_rate  = data_multi_140min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_140min = unique_multifibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_multifibres_140min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_140min_rate_perm4 = len(IDs140min_inBaseline)
perc_multifibre_bundles_140min_rate_perm4 = nr_multifibre_bundles_140min_rate_perm4 / total_nr_multifibre_bundles

#permuation 5
times_140_perm5 = np.arange(80,601,140)
data_multi_140min_rate = data_multi[data_multi['Time'].isin(times_140_perm5)] 
unique_multifibres_140min_rate  = data_multi_140min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_140min = unique_multifibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_multifibres_140min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_140min_rate_perm5 = len(IDs140min_inBaseline)
perc_multifibre_bundles_140min_rate_perm5 = nr_multifibre_bundles_140min_rate_perm5 / total_nr_multifibre_bundles

#permuation 6
times_140_perm6 = np.arange(100,601,140)
data_multi_140min_rate = data_multi[data_multi['Time'].isin(times_140_perm6)] 
unique_multifibres_140min_rate  = data_multi_140min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_140min = unique_multifibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_multifibres_140min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_140min_rate_perm6 = len(IDs140min_inBaseline)
perc_multifibre_bundles_140min_rate_perm6 = nr_multifibre_bundles_140min_rate_perm6 / total_nr_multifibre_bundles

#permuation 7
times_140_perm7 = np.arange(120,601,140)
data_multi_140min_rate = data_multi[data_multi['Time'].isin(times_140_perm7)] 
unique_multifibres_140min_rate  = data_multi_140min_rate.drop_duplicates('TrackID')
list_trackID_multifibres_140min = unique_multifibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_multifibres_140min.isin(list_trackID_multifibres_baseline) 
nr_multifibre_bundles_140min_rate_perm7 = len(IDs140min_inBaseline)
perc_multifibre_bundles_140min_rate_perm7 = nr_multifibre_bundles_140min_rate_perm7 / total_nr_multifibre_bundles

# mean
nr_multifibre_bundles_140min_rate = np.mean([nr_multifibre_bundles_140min_rate_perm1, nr_multifibre_bundles_140min_rate_perm2, 
                                            nr_multifibre_bundles_140min_rate_perm3, nr_multifibre_bundles_140min_rate_perm4,
                                            nr_multifibre_bundles_140min_rate_perm5, nr_multifibre_bundles_140min_rate_perm6,
                                            nr_multifibre_bundles_140min_rate_perm7])
perc_multifibre_bundles_140min_rate = np.mean([perc_multifibre_bundles_140min_rate_perm1, perc_multifibre_bundles_140min_rate_perm2, 
                                              perc_multifibre_bundles_140min_rate_perm3, perc_multifibre_bundles_140min_rate_perm4,
                                              perc_multifibre_bundles_140min_rate_perm5, perc_multifibre_bundles_140min_rate_perm6,
                                              perc_multifibre_bundles_140min_rate_perm7])



frame_hour_rates = [20,40,60,80,100, 120, 140]
nr_detected_multifibre_bundles_perm = [total_nr_multifibre_bundles,
                                  nr_multifibre_bundles_40min_rate, nr_multifibre_bundles_60min_rate,
                                  nr_multifibre_bundles_80min_rate, nr_multifibre_bundles_100min_rate,
                                  nr_multifibre_bundles_120min_rate, nr_multifibre_bundles_140min_rate]
percentage_detected_multifibre_bundles_perm = [perc_total_nr_multifibre_bundles, perc_multifibre_bundles_40min_rate,
                                          perc_multifibre_bundles_60min_rate, perc_multifibre_bundles_80min_rate,
                                          perc_multifibre_bundles_100min_rate, perc_multifibre_bundles_120min_rate,
                                          perc_multifibre_bundles_140min_rate]


data_frame_rates_perm = pd.DataFrame(list(zip(frame_hour_rates, nr_detected_multifibre_bundles_perm, percentage_detected_multifibre_bundles_perm)),
                                columns = ['Frame Rate', 'Nr detected multifibre bundles', 'Percentage detected multifibre bundles'])

x = data_frame_rates_perm['Frame Rate']
y = data_frame_rates_perm['Percentage detected multifibre bundles']
popt, pcov = curve_fit(func, x, y,  p0 = [1, 80, 0.2, 1], maxfev=50000)

# # calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 0 min
# renorm_factor = func(0, *popt)

# calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 5 min
renorm_factor = func(5, *popt)

x_continuous = np.arange(20, 150, 1)
# Plot the frame rate and the number of detected multifibre bundles

x_continuous = np.arange(0, 150, 1)

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
ax.set_xlim([0, 150])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,20,40,60,80,100,120,140])
plt.show(block=True)
 
plt.close()

# Adjusted to renormalization factor

x_continuous = np.arange(0, 150, 1)

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
ax.set_xlim([0, 150])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,20,40,60,80,100,120,140])

#ax.set_xticklabels([0,5,20,40,60,80,100,120,140])
plt.show(block=True)
 
plt.close()


perc_detected_20min = func(20, *popt)/renorm_factor * 100
perc_detected_40min = func(40, *popt)/renorm_factor * 100
difference_perc_detected_20_40 = perc_detected_20min - perc_detected_40min 

print('Percentage of multi-fibre axonal structures being detected with a 40 minute frame rate: ' + str(perc_detected_40min) + '%')
print('Percentage of multi-fibre axonal structures being detected with a 20 minute frame rate: ' + str(perc_detected_20min) + '%')
print('Difference of percentage of multi-fibre axonal structures being detected between 20 and 40 minute frame rate: ' + str(difference_perc_detected_20_40) + '%')



# %% Single-fibre axonal structures

# Times are now from 0 to 600 minutes, with a 20 minute frame rate


########### Considering permutations #######################

# data we get if we measure every 40 minutes

#permutation 1
times_40_perm1 = np.arange(0,601,40)
data_single_40min_rate = data_single[data_single['Time'].isin(times_40_perm1)] 
unique_singlefibres_40min_rate  = data_single_40min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_40min = unique_singlefibres_40min_rate['TrackID']
IDs40min_inBaseline = list_trackID_singlefibres_40min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_40min_rate_perm1 = len(IDs40min_inBaseline)
perc_singlefibre_bundles_40min_rate_perm1 = nr_singlefibre_bundles_40min_rate_perm1 / total_nr_singlefibre_bundles

#permutation 2
times_40_perm2 = np.arange(20,601,40)
data_single_40min_rate = data_single[data_single['Time'].isin(times_40_perm2)] 
unique_singlefibres_40min_rate  = data_single_40min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_40min = unique_singlefibres_40min_rate['TrackID']
IDs40min_inBaseline = list_trackID_singlefibres_40min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_40min_rate_perm2 = len(IDs40min_inBaseline)
perc_singlefibre_bundles_40min_rate_perm2 = nr_singlefibre_bundles_40min_rate_perm2 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_40min_rate = np.mean([nr_singlefibre_bundles_40min_rate_perm1, nr_singlefibre_bundles_40min_rate_perm2])
perc_singlefibre_bundles_40min_rate = np.mean([perc_singlefibre_bundles_40min_rate_perm1, perc_singlefibre_bundles_40min_rate_perm2])


# data we get if we measure every 60 minutes

#permutation 1
times_60_perm1 = np.arange(0,601,60)
data_single_60min_rate = data_single[data_single['Time'].isin(times_60_perm1)] 
unique_singlefibres_60min_rate  = data_single_60min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_60min = unique_singlefibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_singlefibres_60min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_60min_rate_perm1 = len(IDs60min_inBaseline)
perc_singlefibre_bundles_60min_rate_perm1 = nr_singlefibre_bundles_60min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_60_perm2 = np.arange(20,601,60)
data_single_60min_rate = data_single[data_single['Time'].isin(times_60_perm2)] 
unique_singlefibres_60min_rate  = data_single_60min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_60min = unique_singlefibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_singlefibres_60min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_60min_rate_perm2 = len(IDs60min_inBaseline)
perc_singlefibre_bundles_60min_rate_perm2 = nr_singlefibre_bundles_60min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_60_perm3 = np.arange(40,601,60)
data_single_60min_rate = data_single[data_single['Time'].isin(times_60_perm3)] 
unique_singlefibres_60min_rate  = data_single_60min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_60min = unique_singlefibres_60min_rate['TrackID']
IDs60min_inBaseline = list_trackID_singlefibres_60min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_60min_rate_perm3 = len(IDs60min_inBaseline)
perc_singlefibre_bundles_60min_rate_perm3 = nr_singlefibre_bundles_60min_rate_perm3 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_60min_rate = np.mean([nr_singlefibre_bundles_60min_rate_perm1, nr_singlefibre_bundles_60min_rate_perm2, 
                                            nr_singlefibre_bundles_60min_rate_perm3])
perc_singlefibre_bundles_60min_rate = np.mean([perc_singlefibre_bundles_60min_rate_perm1, perc_singlefibre_bundles_60min_rate_perm2, 
                                              perc_singlefibre_bundles_60min_rate_perm3])



# data we get if we measure every 80 minutes

#permutation 1
times_80_perm1 = np.arange(0,601,80)
data_single_80min_rate = data_single[data_single['Time'].isin(times_80_perm1)] 
unique_singlefibres_80min_rate  = data_single_80min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_80min = unique_singlefibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_singlefibres_80min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_80min_rate_perm1 = len(IDs80min_inBaseline)
perc_singlefibre_bundles_80min_rate_perm1 = nr_singlefibre_bundles_80min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_80_perm2 = np.arange(20,601,80)
data_single_80min_rate = data_single[data_single['Time'].isin(times_80_perm2)] 
unique_singlefibres_80min_rate  = data_single_80min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_80min = unique_singlefibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_singlefibres_80min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_80min_rate_perm2 = len(IDs80min_inBaseline)
perc_singlefibre_bundles_80min_rate_perm2 = nr_singlefibre_bundles_80min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_80_perm3 = np.arange(40,601,80)
data_single_80min_rate = data_single[data_single['Time'].isin(times_80_perm3)] 
unique_singlefibres_80min_rate  = data_single_80min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_80min = unique_singlefibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_singlefibres_80min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_80min_rate_perm3 = len(IDs80min_inBaseline)
perc_singlefibre_bundles_80min_rate_perm3 = nr_singlefibre_bundles_80min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_80_perm4 = np.arange(60,601,80)
data_single_80min_rate = data_single[data_single['Time'].isin(times_80_perm4)] 
unique_singlefibres_80min_rate  = data_single_80min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_80min = unique_singlefibres_80min_rate['TrackID']
IDs80min_inBaseline = list_trackID_singlefibres_80min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_80min_rate_perm4 = len(IDs80min_inBaseline)
perc_singlefibre_bundles_80min_rate_perm4 = nr_singlefibre_bundles_80min_rate_perm4 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_80min_rate = np.mean([nr_singlefibre_bundles_80min_rate_perm1, nr_singlefibre_bundles_80min_rate_perm2, 
                                            nr_singlefibre_bundles_80min_rate_perm3, nr_singlefibre_bundles_80min_rate_perm4])
perc_singlefibre_bundles_80min_rate = np.mean([perc_singlefibre_bundles_80min_rate_perm1, perc_singlefibre_bundles_80min_rate_perm2, 
                                              perc_singlefibre_bundles_80min_rate_perm3, perc_singlefibre_bundles_80min_rate_perm4])


# data we get if we measure every 100 minutes

#permutation 1
times_100_perm1 = np.arange(0,601,100)
data_single_100min_rate = data_single[data_single['Time'].isin(times_100_perm1)] 
unique_singlefibres_100min_rate  = data_single_100min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_100min = unique_singlefibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_singlefibres_100min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_100min_rate_perm1 = len(IDs100min_inBaseline)
perc_singlefibre_bundles_100min_rate_perm1 = nr_singlefibre_bundles_100min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_100_perm2 = np.arange(20,601,100)
data_single_100min_rate = data_single[data_single['Time'].isin(times_100_perm2)] 
unique_singlefibres_100min_rate  = data_single_100min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_100min = unique_singlefibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_singlefibres_100min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_100min_rate_perm2 = len(IDs100min_inBaseline)
perc_singlefibre_bundles_100min_rate_perm2 = nr_singlefibre_bundles_100min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_100_perm3 = np.arange(40,601,100)
data_single_100min_rate = data_single[data_single['Time'].isin(times_100_perm3)] 
unique_singlefibres_100min_rate  = data_single_100min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_100min = unique_singlefibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_singlefibres_100min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_100min_rate_perm3 = len(IDs100min_inBaseline)
perc_singlefibre_bundles_100min_rate_perm3 = nr_singlefibre_bundles_100min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_100_perm4 = np.arange(60,601,100)
data_single_100min_rate = data_single[data_single['Time'].isin(times_100_perm4)] 
unique_singlefibres_100min_rate  = data_single_100min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_100min = unique_singlefibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_singlefibres_100min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_100min_rate_perm4 = len(IDs100min_inBaseline)
perc_singlefibre_bundles_100min_rate_perm4 = nr_singlefibre_bundles_100min_rate_perm4 / total_nr_singlefibre_bundles

#permuation 5
times_100_perm5 = np.arange(80,601,100)
data_single_100min_rate = data_single[data_single['Time'].isin(times_100_perm5)] 
unique_singlefibres_100min_rate  = data_single_100min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_100min = unique_singlefibres_100min_rate['TrackID']
IDs100min_inBaseline = list_trackID_singlefibres_100min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_100min_rate_perm5 = len(IDs100min_inBaseline)
perc_singlefibre_bundles_100min_rate_perm5 = nr_singlefibre_bundles_100min_rate_perm5 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_100min_rate = np.mean([nr_singlefibre_bundles_100min_rate_perm1, nr_singlefibre_bundles_100min_rate_perm2, 
                                            nr_singlefibre_bundles_100min_rate_perm3, nr_singlefibre_bundles_100min_rate_perm4,
                                            nr_singlefibre_bundles_100min_rate_perm5])
perc_singlefibre_bundles_100min_rate = np.mean([perc_singlefibre_bundles_100min_rate_perm1, perc_singlefibre_bundles_100min_rate_perm2, 
                                              perc_singlefibre_bundles_100min_rate_perm3, perc_singlefibre_bundles_100min_rate_perm4,
                                              perc_singlefibre_bundles_100min_rate_perm5])

# data we get if we measure every 120 minutes

#permutation 1
times_120_perm1 = np.arange(0,601,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm1)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm1 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm1 = nr_singlefibre_bundles_120min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_120_perm2 = np.arange(20,601,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm2)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm2 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm2 = nr_singlefibre_bundles_120min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_120_perm3 = np.arange(40,601,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm3)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm3 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm3 = nr_singlefibre_bundles_120min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_120_perm4 = np.arange(60,601,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm4)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm4 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm4 = nr_singlefibre_bundles_120min_rate_perm4 / total_nr_singlefibre_bundles

#permuation 5
times_120_perm5 = np.arange(80,601,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm5)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm5 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm5 = nr_singlefibre_bundles_120min_rate_perm5 / total_nr_singlefibre_bundles

#permuation 6
times_120_perm6 = np.arange(100,601,120)
data_single_120min_rate = data_single[data_single['Time'].isin(times_120_perm6)] 
unique_singlefibres_120min_rate  = data_single_120min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_120min = unique_singlefibres_120min_rate['TrackID']
IDs120min_inBaseline = list_trackID_singlefibres_120min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_120min_rate_perm6 = len(IDs120min_inBaseline)
perc_singlefibre_bundles_120min_rate_perm6 = nr_singlefibre_bundles_120min_rate_perm6 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_120min_rate = np.mean([nr_singlefibre_bundles_120min_rate_perm1, nr_singlefibre_bundles_120min_rate_perm2, 
                                            nr_singlefibre_bundles_120min_rate_perm3, nr_singlefibre_bundles_120min_rate_perm4,
                                            nr_singlefibre_bundles_120min_rate_perm5, nr_singlefibre_bundles_120min_rate_perm6])
perc_singlefibre_bundles_120min_rate = np.mean([perc_singlefibre_bundles_120min_rate_perm1, perc_singlefibre_bundles_120min_rate_perm2, 
                                              perc_singlefibre_bundles_120min_rate_perm3, perc_singlefibre_bundles_120min_rate_perm4,
                                              perc_singlefibre_bundles_120min_rate_perm5, perc_singlefibre_bundles_120min_rate_perm6])


# data we get if we measure every 140 minutes

#permutation 1
times_140_perm1 = np.arange(0,601,140)
data_single_140min_rate = data_single[data_single['Time'].isin(times_140_perm1)] 
unique_singlefibres_140min_rate  = data_single_140min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_140min = unique_singlefibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_singlefibres_140min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_140min_rate_perm1 = len(IDs140min_inBaseline)
perc_singlefibre_bundles_140min_rate_perm1 = nr_singlefibre_bundles_140min_rate_perm1 / total_nr_singlefibre_bundles

#permuation 2
times_140_perm2 = np.arange(20,601,140)
data_single_140min_rate = data_single[data_single['Time'].isin(times_140_perm2)] 
unique_singlefibres_140min_rate  = data_single_140min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_140min = unique_singlefibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_singlefibres_140min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_140min_rate_perm2 = len(IDs140min_inBaseline)
perc_singlefibre_bundles_140min_rate_perm2 = nr_singlefibre_bundles_140min_rate_perm2 / total_nr_singlefibre_bundles

#permuation 3
times_140_perm3 = np.arange(40,601,140)
data_single_140min_rate = data_single[data_single['Time'].isin(times_140_perm3)] 
unique_singlefibres_140min_rate  = data_single_140min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_140min = unique_singlefibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_singlefibres_140min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_140min_rate_perm3 = len(IDs140min_inBaseline)
perc_singlefibre_bundles_140min_rate_perm3 = nr_singlefibre_bundles_140min_rate_perm3 / total_nr_singlefibre_bundles

#permuation 4
times_140_perm4 = np.arange(60,601,140)
data_single_140min_rate = data_single[data_single['Time'].isin(times_140_perm4)] 
unique_singlefibres_140min_rate  = data_single_140min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_140min = unique_singlefibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_singlefibres_140min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_140min_rate_perm4 = len(IDs140min_inBaseline)
perc_singlefibre_bundles_140min_rate_perm4 = nr_singlefibre_bundles_140min_rate_perm4 / total_nr_singlefibre_bundles

#permuation 5
times_140_perm5 = np.arange(80,601,140)
data_single_140min_rate = data_single[data_single['Time'].isin(times_140_perm5)] 
unique_singlefibres_140min_rate  = data_single_140min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_140min = unique_singlefibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_singlefibres_140min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_140min_rate_perm5 = len(IDs140min_inBaseline)
perc_singlefibre_bundles_140min_rate_perm5 = nr_singlefibre_bundles_140min_rate_perm5 / total_nr_singlefibre_bundles

#permuation 6
times_140_perm6 = np.arange(100,601,140)
data_single_140min_rate = data_single[data_single['Time'].isin(times_140_perm6)] 
unique_singlefibres_140min_rate  = data_single_140min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_140min = unique_singlefibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_singlefibres_140min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_140min_rate_perm6 = len(IDs140min_inBaseline)
perc_singlefibre_bundles_140min_rate_perm6 = nr_singlefibre_bundles_140min_rate_perm6 / total_nr_singlefibre_bundles

#permuation 7
times_140_perm7 = np.arange(120,601,140)
data_single_140min_rate = data_single[data_single['Time'].isin(times_140_perm7)] 
unique_singlefibres_140min_rate  = data_single_140min_rate.drop_duplicates('TrackID')
list_trackID_singlefibres_140min = unique_singlefibres_140min_rate['TrackID']
IDs140min_inBaseline = list_trackID_singlefibres_140min.isin(list_trackID_singlefibres_baseline) 
nr_singlefibre_bundles_140min_rate_perm7 = len(IDs140min_inBaseline)
perc_singlefibre_bundles_140min_rate_perm7 = nr_singlefibre_bundles_140min_rate_perm7 / total_nr_singlefibre_bundles

# mean
nr_singlefibre_bundles_140min_rate = np.mean([nr_singlefibre_bundles_140min_rate_perm1, nr_singlefibre_bundles_140min_rate_perm2, 
                                            nr_singlefibre_bundles_140min_rate_perm3, nr_singlefibre_bundles_140min_rate_perm4,
                                            nr_singlefibre_bundles_140min_rate_perm5, nr_singlefibre_bundles_140min_rate_perm6,
                                            nr_singlefibre_bundles_140min_rate_perm7])
perc_singlefibre_bundles_140min_rate = np.mean([perc_singlefibre_bundles_140min_rate_perm1, perc_singlefibre_bundles_140min_rate_perm2, 
                                              perc_singlefibre_bundles_140min_rate_perm3, perc_singlefibre_bundles_140min_rate_perm4,
                                              perc_singlefibre_bundles_140min_rate_perm5, perc_singlefibre_bundles_140min_rate_perm6,
                                              perc_singlefibre_bundles_140min_rate_perm7])



frame_hour_rates = [20,40,60,80,100, 120, 140]
nr_detected_singlefibre_bundles_perm = [total_nr_singlefibre_bundles,
                                  nr_singlefibre_bundles_40min_rate, nr_singlefibre_bundles_60min_rate,
                                  nr_singlefibre_bundles_80min_rate, nr_singlefibre_bundles_100min_rate,
                                  nr_singlefibre_bundles_120min_rate, nr_singlefibre_bundles_140min_rate]
percentage_detected_singlefibre_bundles_perm = [perc_total_nr_singlefibre_bundles, perc_singlefibre_bundles_40min_rate,
                                          perc_singlefibre_bundles_60min_rate, perc_singlefibre_bundles_80min_rate,
                                          perc_singlefibre_bundles_100min_rate, perc_singlefibre_bundles_120min_rate,
                                          perc_singlefibre_bundles_140min_rate]


data_frame_rates_perm = pd.DataFrame(list(zip(frame_hour_rates, nr_detected_singlefibre_bundles_perm, percentage_detected_singlefibre_bundles_perm)),
                                columns = ['Frame Rate', 'Nr detected singlefibre bundles', 'Percentage detected singlefibre bundles'])

x = data_frame_rates_perm['Frame Rate']
y = data_frame_rates_perm['Percentage detected singlefibre bundles']
popt, pcov = curve_fit(func, x, y,  p0 = [1, 80, 0.2, 1], maxfev=50000)

# # calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 0 min
# renorm_factor = func(0, *popt)

# calculate renormalization factor so that the proportion of detected fibres is 1 at a frame rate of 5 min
renorm_factor = func(5, *popt)

x_continuous = np.arange(20, 150, 1)
# Plot the frame rate and the number of detected singlefibre bundles

x_continuous = np.arange(0, 150, 1)

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
ax.set_xlim([0, 150])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,20,40,60,80,100,120,140])

plt.show(block=True)
 
plt.close()

# Adjusted to renormalization factor

x_continuous = np.arange(0, 150, 1)

plt.figure(figsize=(9, 7))


fs =18

width = 0.3

plt.rc('grid', linestyle="-", color='lightgrey')
plt.grid(True)
plt.scatter(data_frame_rates_perm['Frame Rate'], data_frame_rates_perm['Percentage detected singlefibre bundles']/renorm_factor, color = 'black', label='Data, renormalized')
plt.plot(x_continuous, func(x_continuous, *popt)/renorm_factor, label='Inverse logistic regression', color='tomato')
plt.legend(fontsize=14)


plt.xlabel('Frame Rate (minutes)', fontsize=fs)
plt.ylabel('Proportion of detected single-fibre axonal structures', fontsize=fs)
#plt.title('Detected singlefibre axonal structures for different frame rates', fontsize = 18)

plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.set_xlim([0, 150])
ax.set_ylim([0, 1.2])
ax.set_xticks([0,5,20,40,60,80,100,120,140])

plt.show(block=True)
 
plt.close()

