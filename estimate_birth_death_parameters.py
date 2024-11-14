import os
import pandas as pd
import copy
import numpy as np
import scipy
from scipy.special import factorial
from scipy.optimize import curve_fit

def kl_divergence(p, q):
    return sum(p[i] *  np.log(p[i]/q[i]) for i in range(len(p)))

def poisson_pmf(k, lambd):
    return (lambd**k/factorial(k)) * np.exp(-lambd)

def poisson_pmf_adjusted(k, lambd):
    return (poisson_pmf(k, lambd) / (1-poisson_pmf(0, lambd)))


def adjust_kl_divergence(k, data_distribution_events, data_distribution_fibres, mean_distribution_fibres):
        
    n_events_list = np.arange(0,7)
    #n_fibers_list = np.arange(1,14)
    n_fibers_list = np.arange(1,9)

    dist_number_fibers = data_distribution_fibres
    mean_nr_fibers = mean_distribution_fibres
    probx = np.zeros(len(n_events_list))
    delta_t = 1
    for i, event_x in enumerate(n_events_list):
        sum_total = 0
        for j, fiber_j in enumerate(n_fibers_list):
            part1 = dist_number_fibers[j]
            part2_nom1 = ((k + fiber_j * k / mean_nr_fibers) * delta_t) ** event_x
            part2_nom2 = np.exp(-((k + fiber_j * k / mean_nr_fibers) * delta_t))
            part2_den = factorial(event_x)
            part2 = part2_nom1 * part2_nom2 / part2_den
            sum_aux = part1 * part2
            sum_total = sum_total + sum_aux
        probx[i] = sum_total
    return (kl_divergence(data_distribution_events, probx))


def fit_function(k, lamb):
    # return poisson.pmf(k, lamb)
    return lamb**k * np.exp(-lamb) / factorial(k) / (1-np.exp(-lamb))

def find_parameters(data_fibres, data_events):

    counts_fibres, bins = np.histogram(data_fibres, bins=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
    
    bins = [1,2,3,4,5,6,7,8]
    
    counts_fibres_density = counts_fibres / np.sum(counts_fibres)
    
    data_distribution_fibres = counts_fibres_density
    data_distribution_fibres[data_distribution_fibres <= 0] = 10 ** -13
    
    parameters, cov_matrix = curve_fit(fit_function, bins, data_distribution_fibres)
    
    counts_events, bins_events = np.histogram(data_events, bins=[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
    counts_events_density = counts_events / np.sum(counts_events)
        
    data_distribution_events = counts_events_density
    data_distribution_events[data_distribution_events <= 0] = 10 ** -13
    
    result_optimization = scipy.optimize.minimize(adjust_kl_divergence, 0.1, args = (data_distribution_events,data_distribution_fibres, parameters[0]))
    k =result_optimization['x'][0]
    
    mean_nr_fibers = parameters[0]
    
    delta = k / mean_nr_fibers
    
    return k, delta, mean_nr_fibers

# Import and process experimental data

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'Fig4E_Fibers-time_TubGFP_Live_updated.xlsx')
data = pd.read_excel(filename, sheet_name = 'Count')


nr_fibres_per_bundle = data['#Fibers'].to_numpy()
stage = data['Stage'].to_numpy()
id_axonal_structure = data['TrackID'].to_numpy()
id_individual =data['ID'].to_numpy()
fate = data['Fate_Posteriori'].to_numpy()
tub_pos = data['#Fiber.TubGFP.pos']

data_df = pd.DataFrame(
    {'ID': id_individual,
     'TrackID': id_axonal_structure,
     'Time': stage,
     'Nr_Fibres': nr_fibres_per_bundle,
     'Fate': fate,
     'N_Tub_Pos' : tub_pos
    })


data_aux = copy.deepcopy(data_df)


data_time_points = [30,31,32,33,34,35,36,37,38,39,40]

# Calculate the numbers of filopodia "events" (one event equals one increase or 
# decrease in the number of filopodia in between two consecutive time points)

# Correct some entries that appear twice
data_sum_bundles = data_aux.groupby(['TrackID','Time', 'Fate']).agg({'Nr_Fibres': 'sum'})
data_sum_bundles = data_sum_bundles.reset_index(level=0)
data_sum_bundles = data_sum_bundles.reset_index(level=0)
data_sum_bundles = data_sum_bundles.reset_index(level=0)
data_aux = data_sum_bundles

diff_list = []
ids_list = data_aux['TrackID'].unique()
for i in ids_list:
    for t in data_time_points:
        subset = data_aux[(data_aux['Time'] == t) & (data_aux['TrackID'] == i)]
        if not subset.empty:
            time_after = t + 1
            subset_next = data_aux[(data_aux['Time'] == time_after) & (data_aux['TrackID'] == i)]
            if not subset_next.empty:
                diff = subset_next.loc[subset_next.index[0], 'Nr_Fibres'] - subset.loc[subset.index[0], 'Nr_Fibres']
                entry_to_add = [i, t, diff]
                diff_list.append(entry_to_add)

diff_df = pd.DataFrame(diff_list, columns = ['TrackID', 'Time', 'N_Fibres'])
diff_df['N_Fibres'] = diff_df['N_Fibres'].abs()
df_events = diff_df


# Numbers of filopodia for time points 30 and 31 together

subset_30 = data_df.loc[data_df['Time'] == 30,:]
subset_31 = data_df.loc[data_df['Time'] == 31,:]
numbers_filo = np.concatenate((np.array(subset_30['Nr_Fibres']), np.array(subset_31['Nr_Fibres'])))

# Number of events for time points 30 and 31 together
 
subset_30_events = df_events.loc[df_events['Time'] == 30,:]
subset_31_events = df_events.loc[df_events['Time'] == 31,:]
numbers_events = np.concatenate((np.array(subset_30_events['N_Fibres']), np.array(subset_31_events['N_Fibres'])))


k, delta, mean_nr_fibres = find_parameters(numbers_filo, numbers_events)

print('Estimated value for k is ' + str(k))
print('Estimated value for delta is ' + str(delta))

