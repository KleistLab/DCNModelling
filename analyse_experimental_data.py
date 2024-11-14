import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.special import factorial
from scipy.optimize import curve_fit
import copy
from scipy.stats import binomtest
import math
plt.style.use('default')

def fit_function(k, lamb):
    return lamb**k * np.exp(-lamb) / factorial(k) / (1-np.exp(-lamb))


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


# Plot the distributions of the number of filopodia per time point along 
# with their Poisson distribution best fit

fs=20
width_bin = 0.95
transperency = 0.7
ksdist = []
times_list = np.arange(30,41,1)
for time_point in times_list:
    subset_time = data_df.loc[data_df['Time'] == time_point,:]
    numbers_filo = np.array(subset_time['Nr_Fibres'])
    bin_edges = np.arange(1,10,1)
    counts, edges = np.histogram(numbers_filo, bins = bin_edges)
    counts_rel = counts/np.sum(counts)
    parameters, cov_matrix = curve_fit(fit_function, bin_edges[0:-1], counts_rel)
    fitted_vals = fit_function(bin_edges[0:-1], *parameters)
    ksdist_val = np.max(np.abs(np.cumsum(fitted_vals)-np.cumsum(counts_rel)))
    ksdist.append(ksdist_val)
    plt.figure(figsize=(9, 7))
    plt.bar(bin_edges[0:-1], counts_rel, width=width_bin, align="center", alpha = transperency,color='grey', label = 'Data')
    plt.plot(
        bin_edges[0:-1],
        fit_function(bin_edges[0:-1], *parameters),
        marker='D', linestyle='-',
        color='k',
        label='Poisson distribution',
        )
    plt.xlabel('Number of filopodia at ' + str(time_point) + ' hAPF', fontsize=fs)
    plt.ylabel('Probability density', fontsize=fs)
    plt.title(str(time_point) + ' hAPF', fontsize=30)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.text(4, 0.25, 'D = ' + str(np.round(ksdist_val,3)), fontsize = 24)
    ax = plt.gca()
    ax.set_xlim([0.5, 8.5])
    ax.set_ylim([0, 0.8])
    plt.legend(fontsize=fs)
    plt.show(block=True)
    
    
    
# Plot with the deviation from poisson fit for each time point
plt.figure(figsize=(9, 7))
plt.plot(times_list, ksdist, marker='o', linestyle='-', color='k')
plt.xlabel('Developmental time (hAPF)', fontsize=fs)
plt.ylabel('Deviation from Poisson fit (D)', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)



# Plot distributions of filopodia "events" and numbers for times 30hAPF and 31hAPF

# Numbers of filopodia for time points 30 and 31 together

subset_30 = data_df.loc[data_df['Time'] == 30,:]
subset_31 = data_df.loc[data_df['Time'] == 31,:]
numbers_filo = np.concatenate((np.array(subset_30['Nr_Fibres']), np.array(subset_31['Nr_Fibres'])))
bin_edges = np.arange(1,10,1)
counts, edges = np.histogram(numbers_filo, bins = bin_edges)
counts_rel = counts/np.sum(counts)
plt.figure(figsize=(9, 7))
plt.bar(bin_edges[0:-1], counts_rel, width=width_bin, align="center", alpha = transperency, color='grey', label = 'Data')
plt.xlabel('Number of filopodia at 30 and 31 hAPF', fontsize=fs)
plt.ylabel('Probability density', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.set_xlim([0.5, 8.5])
ax.set_ylim([0, 0.7])
plt.legend(fontsize=fs)
plt.show(block=True)

# Number of events for time points 30 and 31 together
 
subset_30_events = df_events.loc[df_events['Time'] == 30,:]
subset_31_events = df_events.loc[df_events['Time'] == 31,:]
numbers_events = np.concatenate((np.array(subset_30_events['N_Fibres']), np.array(subset_31_events['N_Fibres'])))
bin_edges = np.arange(0,8,1)
counts, edges = np.histogram(numbers_events, bins = bin_edges)
counts_rel = counts/np.sum(counts)
plt.figure(figsize=(9, 7))
plt.bar(bin_edges[0:-1], counts_rel, width=width_bin, align="center", alpha = transperency,color='grey', label = 'Data')
plt.xlabel('Number of filopodium events at 30 and 31 hAPF', fontsize=fs)
plt.ylabel('Probability density', fontsize=fs)
#plt.title(str(time_point) + ' hAPF', fontsize=30)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
ax.set_xlim([-0.5, 6.5])
ax.set_ylim([0, 0.5])
plt.legend(fontsize=fs)
plt.show(block=True)



# Plot for fibre independence testing + probability of axonal structure
# survival given number of filopodia per structure for data and simulation

data_df_max_fibres = pd.DataFrame(columns=['ID', 'TrackID', 'Time', 'Nr_Fibres', 'Fate'])

id_list = [1]
for i, id_structure in enumerate(data_df['TrackID']):
    if id_structure not in id_list:
        subset = data_df[(data_df['TrackID'] == id_structure)] 
        line_max = subset.loc[subset['Nr_Fibres'].idxmax()]
        data_df_max_fibres = data_df_max_fibres.append(line_max)
        id_list.append(id_structure)
        

max_n_fibers_list = data_df_max_fibres['Nr_Fibres'].unique()
max_n_fibers_list = np.sort(max_n_fibers_list)

prop_selected = np.zeros(len(max_n_fibers_list))

#confidence interval lower bound values
ci_low_75= np.zeros(len(max_n_fibers_list))
#confidence interval higher bound values
ci_high_75 = np.zeros(len(max_n_fibers_list))

for i, n_fibre in enumerate(max_n_fibers_list):
    selected_lines = data_df_max_fibres[(data_df_max_fibres['Nr_Fibres'] == n_fibre) & (data_df_max_fibres['Fate'] == 'SELECTED')]
    transient_lines = data_df_max_fibres[(data_df_max_fibres['Nr_Fibres'] == n_fibre) & (data_df_max_fibres['Fate'] == 'TRANSIENT')]
    n_selected = selected_lines.shape[0]
    n_transient = transient_lines.shape[0]
    n_total = n_selected + n_transient
    prop_selected[i] = n_selected / n_total
    # Calculate confidence intervals with the Clopper–Pearson 
    result = binomtest(k=n_selected, n=n_total)
    ci_interval_75 = np.asarray(result.proportion_ci(method='exact', confidence_level=0.75))
    ci_low_75[i] = ci_interval_75[0]
    ci_high_75[i] = ci_interval_75[1]
    


plt.figure(figsize=(9, 7))
p = np.arange(0.1, 0.91, 0.1)
cmap = plt.get_cmap('Purples', len(p))
for ind, i in enumerate(p):
    plt.plot(np.arange(0, 13.5, 0.001), 1-(1-i)**np.arange(0, 13.5, 0.001), c=cmap(ind))
norm = mpl.colors.Normalize(vmin=0, vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cb = plt.colorbar(sm, ticks=np.linspace(0, 1, 3))
cb.ax.tick_params(labelsize=20)
cb.set_label('$P_{Survival [filopodium]}$', fontsize=fs)
width = 0.3
plt.plot(max_n_fibers_list, prop_selected, 'ko-', c='k', label = 'Data')
plt.legend(loc="best", fontsize=fs)

for i in range(len(max_n_fibers_list)):
    plt.plot([max_n_fibers_list[i] - width/2, max_n_fibers_list[i] - width/2, max_n_fibers_list[i] + width/2, max_n_fibers_list[i] + width/2, max_n_fibers_list[i] - width/2],
              [ci_low_75[i], ci_high_75[i], ci_high_75[i], ci_low_75[i], ci_low_75[i]], c='k')
    # plt.plot([n_fibers_x_axis_data[i] - width/2, n_fibers_x_axis_data[i] + width/2], [(nr_selected/total_nr)[i], (nr_selected/total_nr)[i]], c='r', lw=1)
plt.xlabel('Maximum number of filopodia per axon ', fontsize=fs)
plt.ylabel('$P_{Survival [axon]}$', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close()
