import os
import numpy as np
import plotly.express as px
import math
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from scipy.stats import binomtest
from scipy.special import factorial
from scipy.optimize import curve_fit
plt.style.use('default')

# Note: the word "fibres" is used throughout this
# file with the same meaning as "filopodia". We also
# use the word "bundle" to refer to "axon" or "axonal structure".

# number of fibres, initially seeded
def fit_function(k, lamb):
    return lamb**k * np.exp(-lamb) / factorial(k) / (1-np.exp(-lamb))

def f_inhibition(m_stable_fibers, B50):
    return (B50 / (m_stable_fibers + B50))

# x1 is the number of transient fibers, x2 is the number of stable fibers
def getPropensities (x1, x2, k, delta, gamma, B50):
    r1 = k
    r2 = delta * x1
    r3 =  x1 * gamma  * f_inhibition(x2, B50)

    result = np.array([r1, r2, r3])
    return result 

def sigmoid(x):
    return 1/(1 + np.exp(-(x-4)))


############### Import experimental data ###############

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

# create new column with number of transient fibers (nr of total fibers - number of tubulin positive fibers)
data_df['N_Transient'] = data_df['Nr_Fibres'] - data_df['N_Tub_Pos']

# Fit number the distribution of the number of fibres at times points 30hAPF
# and 31hAPF with a Poisson distribution, and calculate the rate parameter (lambda_poisson) 
# of that fitted Poisson distribution. This will be used to draw the initial
# number of filopodia in the simulation

subset_30 = data_df.loc[data_df['Time'] == 30,:]
subset_31 = data_df.loc[data_df['Time'] == 31,:]
numbers_filo = np.concatenate((np.array(subset_30['Nr_Fibres']), np.array(subset_31['Nr_Fibres'])))
bin_edges = np.arange(1,10,1)
counts, edges = np.histogram(numbers_filo, bins = bin_edges)
counts_rel = counts/np.sum(counts)
parameters_poisson, cov_matrix = curve_fit(fit_function, bin_edges[0:-1], counts_rel)


############### Simulation #################

# Gillespie simulation for found parameters

# Simulation parameters (time of simulation, stabilization rate and B50 parameter
#  of the inhibition funciton)

N_gil = 50000 # number of gillespie simulations
T = 20 # time of simulation
gamma = 0.0135
B50 = 4
lambda_poisson = parameters_poisson[0]
k = 0.520
delta = 0.535
 
# options for the saving of the data
time_disc = 1 # take a snapshot every minute (time discretization)
time_points = np.ceil(T/time_disc)
 
# initialize storage matrix (dimension is nr of simulations x nr of storage time points)
ens_data_x1 = np.zeros((int(N_gil), int(time_points))) # ensemble data transient fibers
ens_data_x2 = np.zeros((int(N_gil), int(time_points))) # ensemble data stable fibers
  
fate_list = np.zeros(N_gil)
max_num_fibres_list = np.zeros(N_gil)
 
 
# start of stochastic simulations
for i in range(0, N_gil):
          
    time_count = 1    
    t = 0
    x1_zero = 0
    max_num_fibres = 0
    while(x1_zero == 0):
        x1_zero = np.random.poisson(lambda_poisson)
    x2_zero = 0
    x1 = x1_zero
    x2 = x2_zero
 
    ens_data_x1[i, 0] = x1
    ens_data_x2[i, 0] = x2
 
    stopped =False
 
     
    #start one growth cone simulation
    while (t < T):
         
        # compute protensities a(X,t)
        a = getPropensities(x1, x2, k, delta, gamma, B50)
        a0 = np.sum(a)
         
        # 1. sample random waiting time tau from exponential distribution
        random1 = np.random.rand()
        #random1 = 0.1
        tau = (1/a0)*math.log(1/random1)  
        t = t+tau # Update time
         
        # 2. sample reaction
        random2 = np.random.rand()
 
        j = np.array(np.where(random2 <= np.cumsum(a)/a0))[0][0]
                  
        #3. Update state according to stoichiometric vector of the choosen
        #reaction
         
        if (time_count <= T):
            if (t > time_count):
                last_t = int(np.floor(t -tau))
                time_count = time_count + 1
         
        if (j == 0):
            x1 = x1 + 1
            x = x1 + x2

 
        elif(j == 1):
            x1 = x1 - 1
            x = x1 + x2

 
        elif(j == 2):
            x2 = x2 + 1
            x1 = x1 - 1
         

        num_fibres_now = x1 + x2
         
        print('time: ' +str(t))
        print('nr_fibres: ' + str(num_fibres_now))
         
        if(num_fibres_now > max_num_fibres):
            max_num_fibres = num_fibres_now

        if(x1 < 0):
            x1 = 0
         
        #---Save simulations at the storage points into ens_Data 
        t_now = int(np.ceil(t))
        t_before = int(np.ceil(t-tau))
        t_passed = np.arange(t_before+1, t_now, time_disc)
        t_passed = t_passed.astype(np.int64)

        if(((x1 == 0) & (x2 == 0))):
            fate = 1 # filopodia is transient
            if(t_now >= T):
                t_now = T-1
                t_passed = np.arange(t_before+1, t_now, time_disc)
                t_passed = t_passed.astype(np.int64)
                 
            ens_data_x1[i, t_now:int(time_points)] = x1
            if (t_passed.size != 0):
               aux_0 = ens_data_x1[i, t_before]
               ens_data_x1[i, t_passed] = aux_0
            ens_data_x2[i, t_now:int(time_points)] = x2
            if (t_passed.size != 0):
               aux_1 = ens_data_x2[i, t_before]
               ens_data_x2[i, t_passed] = aux_1
            stopped = True
            break
         
        if (x2 > 0):
            fate = 2 # filopodia is now stable
             
        if (t_now >= T):
            break
         
        #update ensemble matrices (storage)        
        ens_data_x1[i, t_now] = x1
        ens_data_x2[i, t_now] = x2
        
        # if t_passed isn't 0, fill in data for the time between
        if (t_passed.size != 0):
            
           aux_0 = ens_data_x1[i, t_before]
           ens_data_x1[i, t_passed] = aux_0
            
           aux_1 = ens_data_x2[i, t_before]
           ens_data_x2[i, t_passed] = aux_1
             
    fate_list[i] = fate
  
    max_num_fibres_list[i] = max_num_fibres
    if (stopped == False)   :   
        #fill in last bits of ensemble matrices
        t_before = int(np.ceil(t-tau))
        t_passed = np.arange(t_before+1, T, time_disc)
        t_passed = t_passed.astype(np.int64)
         
        if (t_passed.size != 0):
             
            aux_0 = ens_data_x1[i, t_before]
            ens_data_x1[i, t_passed] = aux_0
            aux_1 = ens_data_x2[i, t_before]
            ens_data_x2[i, t_passed] = aux_1
             

      
         
############### ANALYSE SIMULATION RESULTS #######################

zero_positions = np.where(max_num_fibres_list == 0)[0]
fate_list = np.delete(fate_list, zero_positions)
max_num_fibres_list = np.delete(max_num_fibres_list, zero_positions)
 

# ensemble data for x- count index of first entry that is zero. those zero entries are deleted
df_fibres = pd.DataFrame(columns=['Time','Nr_Fibres', 'ID'])

ens_data_x = ens_data_x1 + ens_data_x2
 
for i in range(len(ens_data_x)):
    if 0 in ens_data_x[i]:
        id_first_zero_x = int(np.argwhere(ens_data_x[i] == 0)[0])
        x_without_zero = ens_data_x[i,0:id_first_zero_x]
    else: 
        x_without_zero = ens_data_x[i,:]
     
    time_x = np.arange(30, 30+len(x_without_zero), 1)
    id_list_x = [i] * len(x_without_zero)
    new_df_x = pd.DataFrame(zip(time_x, x_without_zero, id_list_x), columns=['Time','Nr_Fibres', 'ID'])
    df_fibres = df_fibres.append(new_df_x, ignore_index=True)
 
 
last_time_point = ens_data_x2[:,-1]
 
selected = last_time_point > 0
 
total = len(selected)
count_selected = sum(selected)
count_not_selected = total - count_selected
perc_selected = count_selected / total 
perc_not_selected = count_not_selected / total
 
# different way of calculaitng the same thing:
perc_selected_2 = sum(fate_list==2)/N_gil * 100
 
# consider now only the ones that were selected
 
only_selected = last_time_point[last_time_point > 0]
 
unique, counts = np.unique(only_selected, return_counts=True)
  
sum_counts = len(only_selected)
 
counts_simulation = counts/sum_counts 
counts_proportion = counts/sum_counts * 100
 
x_vals= ['Selected', 'Transient']
y_vals = [perc_selected, perc_not_selected]
labels_perc = []
for i in range(len(y_vals)):
    perc = y_vals[i] * 100
    text_to_add = str(np.round(perc)) + "%" 
    labels_perc.append(text_to_add)
 
# Bar plot for percentage of transient vs stable (medulla targetting) axonal structures
plt.figure(figsize=(9, 7))
fs =20
width = 0.3
y_vals_perc = y_vals * 100
bars = plt.bar(x_vals, y_vals, color='dimgray')
A_as_ticklabel = [f"{100*a:.1f}%" for a in y_vals]
plt.bar_label(bars, labels=A_as_ticklabel, size = fs)
plt.xlabel('Fate of axon', fontsize=fs)
plt.ylabel('Proportion of M-DCN axons', fontsize=fs)
plt.suptitle('Simulation: Proportion of stable vs transient axons ', fontsize = fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close()

x_vals= unique
y_vals = counts_simulation
labels_perc = []
for i in range(len(y_vals)):
    perc = y_vals[i] * 100
    text_to_add = str(np.round(perc)) + "%" 
    labels_perc.append(text_to_add)

# Remove bar if bar smaller than 0.5%
small_values = np.where(y_vals < 0.005)
x_vals = np.delete(x_vals, small_values)
y_vals= np.delete(y_vals, small_values)


# Bar plot for number of filopodia per medulla targetting axonal structure 
plt.figure(figsize=(9, 7))
fs =20
width = 0.3
bars = plt.bar(x_vals, y_vals, color='dimgray')
A_as_ticklabel = [f"{100*a:.00f}%" for a in y_vals]
plt.bar_label(bars, labels=A_as_ticklabel, size = fs)
plt.xlabel('Number of filopodia', fontsize=fs)
plt.ylabel('Proportion of M-DCN axons', fontsize=fs)
plt.suptitle('Simulation: Number of filopodia per stable (M-DCN) axon',  fontsize = fs)
plt.xticks(x_vals,fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close()
 
# Convert simulation results to a same timescale
 
tids = np.array(data_df['TrackID'])
time = np.array(data_df['Time'])
first_time = []
for tid in set(tids):
    currtime_all = time[tid == tids]
    first_time.append(currtime_all[0])
     
bins_histo = [29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5, 38.5, 39.5, 40.5]
counts, bins = np.histogram(first_time, bins=bins_histo)
counts_density = counts / np.sum(counts)
 
bins_plot = [30,31,32,33,34,35,36,37,38,39,40]
  
 
cumulative_dist = np.cumsum(counts_density)
 
 
N = N_gil
 
sampled_times = []
for i in range(N):
 
    random2 = np.random.rand()
 
    j = np.array(np.where(random2 <= cumulative_dist))[0][0]
 
    sampled_value = bins_plot[j]
          
    sampled_times.append(sampled_value)
 
 
 
counts, bins = np.histogram(sampled_times, bins=bins_histo)
counts_density = counts / np.sum(counts)
 
bins_plot = [30,31,32,33,34,35,36,37,38,39,40]

 
for id_structure in range(N_gil):
    #subset df for the structure with such id
    #correct times for that subset
    df_fibres.loc[df_fibres["ID"]==id_structure, "Time"] = df_fibres.loc[df_fibres["ID"]==id_structure, "Time"] + sampled_times[id_structure] -30
 
 
 
# Plot the probability of axon survival given maximum number of filopodia per axon

# Data

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
    
 
# percentage of bundles selected and not selected
 
fate_column = data_df[["TrackID", "Fate"]]
fate_column_distinct = fate_column.drop_duplicates()

# Correct case for TrackID 1100000011. This bundle should have fate 'TRANSIENT'
#fate_column_distinct = fate_column_distinct[(fate_column_distinct['Fate'] != 'SELECTED') & (fate_column_distinct['TrackID'] != 1100000011)]
fate_column_distinct = fate_column_distinct.drop(fate_column_distinct[(fate_column_distinct.Fate == 'SELECTED') & (fate_column_distinct.TrackID == 1100000011)].index)
 
 
# Simulation

# Calculate for each bundle in the simulation, the maximum number of fibres
# that bundle has throughout time, and the fate of the bundle
df_fibres_fate = pd.DataFrame({'Max_Nr_Fibres':max_num_fibres_list, 'Fate_Bundle':fate_list})
 
max_n_fibers_list_sim = df_fibres_fate['Max_Nr_Fibres'].unique()
max_n_fibers_list_sim = np.sort(max_n_fibers_list_sim)
 
prop_selected_sim = np.zeros(len(max_n_fibers_list_sim))
 
prop_selected_sim = []
for n_fibre in max_n_fibers_list_sim:
    selected_lines = df_fibres_fate[(df_fibres_fate['Max_Nr_Fibres'] == n_fibre) & (df_fibres_fate['Fate_Bundle'] == 2)]
    transient_lines = df_fibres_fate[(df_fibres_fate['Max_Nr_Fibres'] == n_fibre) & (df_fibres_fate['Fate_Bundle'] == 1)]

    n_selected = selected_lines.shape[0]
    n_transient = transient_lines.shape[0]
    n_total = n_selected + n_transient
    prop_selected_sim.append(n_selected / n_total)
     
 
n_fibers_axis_sim = max_n_fibers_list_sim
n_fibers_axis_theoretical = np.arange(0, max_n_fibers_list_sim[-1]+1,1)
 

# Plot simulation and data curves together
 
plt.figure(figsize=(9, 7))
p = np.arange(0.1, 0.91, 0.1)
cmap = plt.get_cmap('Purples', len(p))
norm = mpl.colors.Normalize(vmin=0, vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
width = 0.3
plt.plot(max_n_fibers_list_sim, prop_selected_sim, 'ko-', c ='r', label = 'Simulation')
plt.plot(max_n_fibers_list, prop_selected, 'ko-', c='k', label = 'Data')
plt.legend(loc="lower right", fontsize=fs)
for i in range(len(max_n_fibers_list)):
    plt.plot([max_n_fibers_list[i] - width/2, max_n_fibers_list[i] - width/2, max_n_fibers_list[i] + width/2, max_n_fibers_list[i] + width/2, max_n_fibers_list[i] - width/2],
              [ci_low_75[i], ci_high_75[i], ci_high_75[i], ci_low_75[i], ci_low_75[i]], c='k')
plt.xlabel('Maximum number of filopodia per axon', fontsize=fs)
plt.ylabel('$P_{Survival [axon]}$', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True) 
plt.close()
 
   
# Plot simulation, data and filopodia independence testing curves together

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
plt.plot(max_n_fibers_list_sim, prop_selected_sim, 'ko-', c ='r', label = 'Simulation')
plt.legend(loc="lower right", fontsize=fs)

for i in range(len(max_n_fibers_list)):
    plt.plot([max_n_fibers_list[i] - width/2, max_n_fibers_list[i] - width/2, max_n_fibers_list[i] + width/2, max_n_fibers_list[i] + width/2, max_n_fibers_list[i] - width/2],
              [ci_low_75[i], ci_high_75[i], ci_high_75[i], ci_low_75[i], ci_low_75[i]], c='k')
plt.xlabel('Maximum number of filopodia per axon ', fontsize=fs)
plt.ylabel('$P_{Survival [axon]}$', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close()


 