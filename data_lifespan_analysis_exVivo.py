import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('default')


# Import data
dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'dataset-20231211.csv')
data = pd.read_csv(filename, sep = ';')
data = data.iloc[:, :-1]

# Remove reference lines from the data
data = data[data.Name == 'structure']

# replace commas by dots in column "Stage"
data["Stage"]=data["Stage"].str.replace(',','.')


# Keep only multi-fibre structures
data_multi = data[data.Position_Type_Track == 'multi']

# Keep only single-fibre structures
data_single = data[data.Position_Type_Track == 'single']

# Calculate initial time, last time, and lifetime of each axonal structure

# In hours
data_multi_lifetime_hours = data_multi.groupby("TrackID")['Stage'].agg(['first', 'last'])
data_multi_lifetime_hours["first"] = pd.to_numeric(data_multi_lifetime_hours["first"])
data_multi_lifetime_hours["last"] = pd.to_numeric(data_multi_lifetime_hours["last"])
data_multi_lifetime_hours["lifetime"] = data_multi_lifetime_hours["last"] - data_multi_lifetime_hours["first"]

data_single_lifetime_hours = data_single.groupby("TrackID")['Stage'].agg(['first', 'last'])
data_single_lifetime_hours["first"] = pd.to_numeric(data_single_lifetime_hours["first"])
data_single_lifetime_hours["last"] = pd.to_numeric(data_single_lifetime_hours["last"])
data_single_lifetime_hours["lifetime"] = data_single_lifetime_hours["last"] - data_single_lifetime_hours["first"]


# In minutes
data_multi_lifetime_minutes = data_multi.groupby("TrackID")['Time'].agg(['first', 'last'])
data_multi_lifetime_minutes["first"] = pd.to_numeric(data_multi_lifetime_minutes["first"])
data_multi_lifetime_minutes["last"] = pd.to_numeric(data_multi_lifetime_minutes["last"])
data_multi_lifetime_minutes["lifetime"] = data_multi_lifetime_minutes["last"] - data_multi_lifetime_minutes["first"]

data_single_lifetime_minutes = data_single.groupby("TrackID")['Time'].agg(['first', 'last'])
data_single_lifetime_minutes["first"] = pd.to_numeric(data_single_lifetime_minutes["first"])
data_single_lifetime_minutes["last"] = pd.to_numeric(data_single_lifetime_minutes["last"])
data_single_lifetime_minutes["lifetime"] = data_single_lifetime_minutes["last"] - data_single_lifetime_minutes["first"]


# Histograms of lifetime(hours)

# Multi-fibre
plt.figure(figsize=(9, 7))
fs =18
width = 1
hist, edges = np.histogram(data_multi_lifetime_hours['lifetime'], bins = np.arange(0,11,1))
freq = hist / float(hist.sum())
plt.bar(edges[0:-1], freq, width=width, align="edge", color='darkgrey')
#plt.legend(fontsize=14)
plt.xlabel('Lifetime of axonal structures (hours)', fontsize=fs)
plt.ylabel('Probability density', fontsize=fs)
plt.title('Multi-fibre axonal structures', fontsize = 18)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
# ax.set_xlim([0, 150])
ax.set_ylim([0, 1])
ax.set_xticks(np.arange(0,11))
plt.show(block=True)


# Single-fibre
plt.figure(figsize=(9, 7))
fs =18
width = 1
hist, edges = np.histogram(data_single_lifetime_hours['lifetime'], bins = np.arange(0,11,1))
freq = hist / float(hist.sum())
plt.bar(edges[0:-1], freq, width=width, align="edge", color='darkgrey')
#plt.legend(fontsize=14)
plt.xlabel('Lifetime of axonal structures (hours)', fontsize=fs)
plt.ylabel('Probability density', fontsize=fs)
plt.title('Single-fibre axonal structures', fontsize = 18)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
ax = plt.gca()
# ax.set_xlim([0, 150])
ax.set_ylim([0, 1])
ax.set_xticks(np.arange(0,11))
plt.show(block=True)

# Histograms of lifetime(minutes)

# Multi-fibre
plt.figure(figsize=(9, 7))
fs =18
width = 20
hist, edges = np.histogram(data_multi_lifetime_minutes['lifetime'], bins = np.arange(0,601,20))
freq = hist / float(hist.sum())
plt.bar(edges[0:-1], freq, width=width, align="edge", color='darkgrey')
#plt.legend(fontsize=14)
plt.xlabel('Lifetime of axonal structures (minutes)', fontsize=fs)
plt.ylabel('Probability density', fontsize=fs)
plt.title('Multi-fibre axonal structures', fontsize = 18)
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
# ax.set_xlim([0, 150])
ax.set_ylim([0, 1])
ax.set_xticks(np.arange(0,601,20))
ax.tick_params(axis='x', labelrotation=90)
plt.show(block=True)

# Single-fibre
plt.figure(figsize=(9, 7))
fs =18
width = 20
hist, edges = np.histogram(data_single_lifetime_minutes['lifetime'], bins = np.arange(0,601,20))
freq = hist / float(hist.sum())
plt.bar(edges[0:-1], freq, width=width, align="edge", color='darkgrey')
#plt.legend(fontsize=14)
plt.xlabel('Lifetime of axonal structures (minutes)', fontsize=fs)
plt.ylabel('Probability density', fontsize=fs)
plt.title('Single-fibre axonal structures', fontsize = 18)
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
# ax.set_xlim([0, 150])
ax.set_ylim([0, 1])
ax.set_xticks(np.arange(0,601,20))
ax.tick_params(axis='x', labelrotation=90)
plt.show(block=True)


# Find unique trackID values for both single and multi-fibre structures

data_multi_lifetime_hours = data_multi_lifetime_hours.reset_index()
data_multi_lifetime_minutes = data_multi_lifetime_minutes.reset_index()
data_single_lifetime_hours = data_single_lifetime_hours.reset_index()
data_single_lifetime_minutes = data_single_lifetime_minutes.reset_index()

data_multi_unique_hours  = data_multi_lifetime_hours.drop_duplicates('TrackID')
data_multi_unique_minutes  = data_multi_lifetime_minutes.drop_duplicates('TrackID')

data_single_unique_hours  = data_single_lifetime_hours.drop_duplicates('TrackID')
data_single_unique_minutes  = data_single_lifetime_minutes.drop_duplicates('TrackID')

# Assign new IDs from 1 to max number of axonal structures
data_multi_unique_hours['ID'] = np.arange(1, len(data_multi_unique_hours)+1,1)
data_multi_unique_minutes['ID'] = np.arange(1, len(data_multi_unique_minutes)+1,1)
data_single_unique_hours['ID'] = np.arange(1, len(data_single_unique_hours)+1,1)
data_single_unique_minutes['ID'] = np.arange(1, len(data_single_unique_minutes)+1,1)


# Plot trajectories in minutes

# Multi-fibres
plt.figure(figsize=(9, 7))
fs =18
width = 60
for i in range(len(data_multi_unique_hours)):
    x1 = data_multi_unique_hours.loc[data_multi_unique_hours.index[i], 'first']
    x2 = data_multi_unique_hours.loc[data_multi_unique_hours.index[i], 'last']
    y1 = data_multi_unique_hours.loc[data_multi_unique_hours.index[i], 'ID']
    y2 = y1
    if (x1 == x2):
        plt.scatter(x1,y1, color='steelblue', marker='.', s=20)
    else:
        plt.plot((x1, x2), (y1, y2), color='steelblue', linewidth=2)
#plt.legend(fontsize=14)
plt.xlabel('Developmental time (hours)', fontsize=fs)
plt.ylabel('Axonal structure ID', fontsize=fs)
plt.title('Lifespan trajectory of multi-fibre axonal structures', fontsize = 18, y = 1.05 )
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.text(30,105.7, 'n=' + str(len(data_multi_unique_hours)), fontsize=fs)
ax = plt.gca()
ax.set_xlim([30, 40])
ax.set_ylim([0, len(data_multi_unique_hours)+1])
#ax.set_xticks([0,5,20,40,60,80,100,120,140])
plt.show(block=True)
 
plt.close()

# Single-fibres
plt.figure(figsize=(9, 7))
fs =18
width = 60
for i in range(len(data_single_unique_hours)):
    x1 = data_single_unique_hours.loc[data_single_unique_hours.index[i], 'first']
    x2 = data_single_unique_hours.loc[data_single_unique_hours.index[i], 'last']
    y1 = data_single_unique_hours.loc[data_single_unique_hours.index[i], 'ID']
    y2 = y1
    if (x1 == x2):
        plt.scatter(x1,y1, color='steelblue', marker='.', s=20)
    else:
        plt.plot((x1, x2), (y1, y2), color='steelblue', linewidth=2)
#plt.legend(fontsize=14)
plt.xlabel('Developmental time (hours)', fontsize=fs)
plt.ylabel('Axonal structure ID', fontsize=fs)
plt.title('Lifespan trajectory of single-fibre axonal structures', fontsize = 18, y=1.05)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.text(30,88.7, 'n=' + str(len(data_single_unique_hours)), fontsize=fs)
ax = plt.gca()
ax.set_xlim([30, 40])
ax.set_ylim([0, len(data_single_unique_hours)+1])
#ax.set_xticks([0,5,20,40,60,80,100,120,140])
plt.show(block=True)
 
plt.close()


# Calculate mean and median lifetimes of axonal structures

# Mean in hours
mean_multifibre_hours = np.mean(data_multi_unique_hours['lifetime'])
mean_singlefibre_hours = np.mean(data_single_unique_hours['lifetime'])


x_vals = ['Multifibre', 'Single Fibre']
y_vals = [mean_multifibre_hours, mean_singlefibre_hours]

plt.figure(figsize=(9, 7))
fs =20
width = 0.3
bars = plt.bar(x_vals, y_vals, color='dimgray')
A_as_ticklabel = [f"{a:.1f}" for a in y_vals]
plt.bar_label(bars, labels=A_as_ticklabel, size = fs)
plt.xlabel('Type of axonal structure', fontsize=fs)
plt.ylabel('Mean  lifetime (hours)', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close()

# Mean in minutes
mean_multifibre_minutes = np.mean(data_multi_unique_minutes['lifetime'])
mean_singlefibre_minutes = np.mean(data_single_unique_minutes['lifetime'])


x_vals = ['Multifibre', 'Single Fibre']
y_vals = [mean_multifibre_minutes, mean_singlefibre_minutes]

plt.figure(figsize=(9, 7))
fs =20
width = 0.3
bars = plt.bar(x_vals, y_vals, color='dimgray')
A_as_ticklabel = [f"{a:.1f}" for a in y_vals]
plt.bar_label(bars, labels=A_as_ticklabel, size = fs)
plt.xlabel('Type of axonal structure', fontsize=fs)
plt.ylabel('Mean  lifetime (minutes)', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close()



# Median in hours
median_multifibre_hours = np.median(data_multi_unique_hours['lifetime'])
median_singlefibre_hours = np.median(data_single_unique_hours['lifetime'])


x_vals = ['Multi-fibre', 'Single-fibre']
y_vals = [median_multifibre_hours, median_singlefibre_hours]

plt.figure(figsize=(9, 7))
fs =20
width = 0.3
bars = plt.bar(x_vals, y_vals, color='dimgray')
A_as_ticklabel = [f"{a:.1f}" for a in y_vals]
plt.bar_label(bars, labels=A_as_ticklabel, size = fs)
plt.xlabel('Type of axonal structure', fontsize=fs)
plt.ylabel('Median  lifetime (hours)', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close() 


# Median in minutes
median_multifibre_minutes = np.median(data_multi_unique_minutes['lifetime'])
median_singlefibre_minutes = np.median(data_single_unique_minutes['lifetime'])


x_vals = ['Multi-fibre', 'Single-fibre']
y_vals = [median_multifibre_minutes, median_singlefibre_minutes]

plt.figure(figsize=(9, 7))
fs =20
width = 0.3
bars = plt.bar(x_vals, y_vals, color='dimgray')
A_as_ticklabel = [f"{a:.1f}" for a in y_vals]
plt.bar_label(bars, labels=A_as_ticklabel, size = fs)
plt.xlabel('Type of axonal structure', fontsize=fs)
plt.ylabel('Median  lifetime (minutes)', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show(block=True)
plt.close() 


# Calculate percentage of multi-fibre axonal structures that have a lifespan of 20 minutes or more:

nr_multi_lifetime_more_equal_20 = data_multi_unique_minutes[data_multi_unique_minutes.lifetime >= 20].shape[0]
perc_multi_lifetime_more_equal_20 = nr_multi_lifetime_more_equal_20 / len(data_multi_unique_minutes) * 100

print('Percentage of multi-fibre axonal structures that have a lifespan of 20 minutes or more (>=20): ' + str(perc_multi_lifetime_more_equal_20) + '%')

# Calculate percentage of multi-fibre axonal structures that have of more than 20 minutes:

nr_multi_lifetime_more20 = data_multi_unique_minutes[data_multi_unique_minutes.lifetime > 20].shape[0]
perc_multi_lifetime_more20 = nr_multi_lifetime_more20 / len(data_multi_unique_minutes) * 100

print('Percentage of multi-fibre axonal structures that have a lifespan of more than 20 minutes (>20): ' + str(perc_multi_lifetime_more20) + '%')


# Calculate percentage of single-fibre axonal structures that have a lifespan of 20 minutes or less:

nr_single_lifetime_less_equal_20 = data_single_unique_minutes[data_single_unique_minutes.lifetime <= 20].shape[0]
perc_single_lifetime_less_equal_20 = nr_single_lifetime_less_equal_20 / len(data_single_unique_minutes) * 100

print('Percentage of single-fibre axonal structures that have a lifespan of 20 minutes or less (<=20): ' + str(perc_single_lifetime_less_equal_20) + '%')


# Calculate percentage of single-fibre axonal structures that have a lifespan of less than 20 minutes:

nr_single_lifetime_less20 = data_single_unique_minutes[data_single_unique_minutes.lifetime < 20].shape[0]
perc_single_lifetime_less20 = nr_single_lifetime_less20 / len(data_single_unique_minutes) * 100

print('Percentage of single-fibre axonal structures that have a lifespan of 20 minutes or less (<20): ' + str(perc_single_lifetime_less20) + '%')


# Calculate percentage of multi-fibre axonal structures that have a lifespan of 30 minutes or more:

nr_multi_lifetime_more_equal_30 = data_multi_unique_minutes[data_multi_unique_minutes.lifetime >= 30].shape[0]
perc_multi_lifetime_more_equal_30 = nr_multi_lifetime_more_equal_30 / len(data_multi_unique_minutes) * 100

print('Percentage of multi-fibre axonal structures that have a lifespan of 30 minutes or more (>=30): ' + str(perc_multi_lifetime_more_equal_30) + '%')

# Calculate percentage of multi-fibre axonal structures that have of more than 30 minutes:

nr_multi_lifetime_more30 = data_multi_unique_minutes[data_multi_unique_minutes.lifetime > 30].shape[0]
perc_multi_lifetime_more30 = nr_multi_lifetime_more30 / len(data_multi_unique_minutes) * 100

print('Percentage of multi-fibre axonal structures that have a lifespan of more than 30 minutes (>30): ' + str(perc_multi_lifetime_more30) + '%')


# Calculate percentage of single-fibre axonal structures that have a lifespan of 30 minutes or less:

nr_single_lifetime_less_equal_30 = data_single_unique_minutes[data_single_unique_minutes.lifetime <= 30].shape[0]
perc_single_lifetime_less_equal_30 = nr_single_lifetime_less_equal_30 / len(data_single_unique_minutes) * 100

print('Percentage of single-fibre axonal structures that have a lifespan of 30 minutes or less (<=30): ' + str(perc_single_lifetime_less_equal_30) + '%')


# Calculate percentage of single-fibre axonal structures that have a lifespan of less than 30 minutes:

nr_single_lifetime_less30 = data_single_unique_minutes[data_single_unique_minutes.lifetime < 30].shape[0]
perc_single_lifetime_less30 = nr_single_lifetime_less30 / len(data_single_unique_minutes) * 100

print('Percentage of single-fibre axonal structures that have a lifespan of 30 minutes or less (<30): ' + str(perc_single_lifetime_less30) + '%')

