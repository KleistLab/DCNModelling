#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 11:10:37 2023

@author: eric
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from statsmodels.stats.proportion import proportions_ztest
#mpl.use('Qt5Agg') # turn on GUI for plotting
import os

dirname = os.path.dirname(__file__)
filename_dat_old = os.path.join(dirname, 'dataset-20231211.csv')
#filename_dat_1 = os.path.join(dirname, 'data_TubGFP-TOM_Culture.xlsx')


# first version, before revision, one animal
dat_old = pd.read_csv(filename_dat_old, sep=';', engine='python')

# second version, after revision, two animals
#dat_1 = pd.read_excel(filename_dat_1, sheet_name='20190725_Filaments')
# dat_2 = pd.read_excel('/home/eric/fibre_stabilization/revision_analysis/data_TubGFP-TOM_Culture.xlsx', sheet_name='20190724_Filaments')
# Maheva: dat_2 and dat_old are the same

# five spatial bins
locs = ['D-Dorsal', 'Dorsal', 'Medial', 'Ventral', 'V-Ventral']

ids_old = dat_old['TrackID']
#ids_1 = dat_1['TrackID']
# ids_2 = dat_2['TrackID']

# get a list of unique IDs
names = dat_old['Name']
ids_uniq_old = np.array(list(set(ids_old[names == 'structure']))) # maybe remove 'F01'
ids_uniq_old = ids_uniq_old[ids_uniq_old != 'F01']

#ids_uniq_1 = np.array(list(set(ids_1)))
# ids_uniq_2 = np.array(list(set(ids_2)))
# ids_uniq_2 = ids_uniq_2[ids_uniq_2 != 'F01']

# loop through the IDs and collect booleans for Tub+ and LP-reaching
LP_bool_old = np.zeros(len(ids_uniq_old), dtype='bool')
tub_bool_old = np.zeros(len(ids_uniq_old), dtype='bool')
reg = []
for idd in ids_uniq_old:
    print('\n' + str(idd))
    print(np.array(dat_old['Location'][ids_old == idd]))
    # in some cases, the Location column contains mixed locations for a given TrackID
    reg.append(np.array(dat_old['Location'][ids_old == idd])[0])
    for ll in np.array(dat_old['Length'][ids_old == idd]):
        print(ll)
        if type(ll)==str:
            if 'LP' in ll:
                LP_bool_old[ids_uniq_old==idd] = True
    for tt in np.array(dat_old['Position.TubGFP'][ids_old == idd]):
        print(tt)
        if type(tt)==str:
            if 'pos' in tt:
                tub_bool_old[ids_uniq_old==idd] = True

# # same thing for the second data set
# LP_bool_1 = np.zeros(len(ids_uniq_1), dtype='bool')
# tub_bool_1 = np.zeros(len(ids_uniq_1), dtype='bool')
# for idd in ids_uniq_1:
#     print('\n' + str(idd))
#     print(np.array(dat_1['Location'][ids_1 == idd]))
#     # in some cases, the Location column contains mixed locations for a given TrackID
#     reg.append(np.array(dat_1['Location'][ids_1 == idd])[0])
#     for ll in np.array(dat_1['Length'][ids_1 == idd]):
#         print(ll)
#         if type(ll)==str:
#             if 'LP' in ll:
#                 LP_bool_1[ids_uniq_1==idd] = True
#     for tt in np.array(dat_1['Position.TubGFP'][ids_1 == idd]):
#         print(tt)
#         if type(tt)==str:
#             if 'pos' in tt:
#                 tub_bool_1[ids_uniq_1==idd] = True

# LP_bool = np.r_[LP_bool_old, LP_bool_1]
# tub_bool = np.r_[tub_bool_old, tub_bool_1]
reg = np.array(reg)

LP_bool = LP_bool_old
tub_bool = tub_bool_old

# for lobula-plate reaching axons
# conf_low = np.array([0.0192, 0.0259, 0.0192, 0.0259, 0.0016])
# conf_high = np.array([0.4545, 0.1730, 0.2433, 0.1730, 0.3023])
conf_low = np.array([0.3686, 0.2083, 0.1015, 0.1313, 0.2444]) # Pearson-Clopper confidence intervals, calculated in an external spreadsheet
conf_high = np.array([0.5638, 0.3752, 0.2788, 0.2774, 0.4537])

# calculate proportions and plot results
counts = []
nobs = []
plt.figure()
for i, loc in enumerate(locs):
    mea = np.mean(LP_bool[reg==loc])
    plt.bar(i, mea, ec='k', lw=2, fc='None')
    #plt.errorbar(i, mea, yerr=[[mea-conf_low[i]], [conf_high[i]-mea]], color='k')
    count = int(np.sum(LP_bool[reg==loc]))
    counts.append(count)
    nob = int(np.sum(reg==loc))
    nobs.append(nob)
    plt.text(i-0.04, np.mean(LP_bool[reg==loc])+0.005, str(count)+'/'+str(nob), horizontalalignment='right')
    plt.text(i+0.1, np.mean(LP_bool[reg==loc])+0.005, 'n.s.', horizontalalignment='left')
plt.xticks(range(5), locs)
plt.ylabel('Proportion of LP reaching axons')
plt.savefig(str(dirname)+'/LP_reaching_axons.png')

print('one sided')
for i in range(5):
    print(np.sum(counts[:i] + counts[i+1 :])/np.sum(nobs[:i] + nobs[i+1 :]))
    print(proportions_ztest(counts[i], nobs[i], value=np.sum(counts[:i] + counts[i+1 :])/np.sum(nobs[:i] + nobs[i+1 :])))

print('two sided')
for i in range(5):
    for j in range(5):
        if j > i:
            print(proportions_ztest([counts[i], counts[j]], [nobs[i], nobs[j]]))


# for tubulin positive fibers
# conf_low = np.array([0.0019, 0.0098, 0.034, 0.0259, 0.016]) # dat_old
# conf_high = np.array([0.3603, 0.1309, 0.282, 0.173, 0.3023])

# conf_low = np.array([0.0259, 0.0227, 0.0007, 0.0123, 0.032]) # dat_1
# conf_high = np.array([0.2214, 0.1960, 0.1533, 0.1624, 0.2674])

conf_low = np.array([0.0296, 0.0253, 0.0247, 0.0305, 0.0326]) # dat_old + dat_1
conf_high = np.array([0.1962, 0.1235, 0.1656, 0.1325, 0.2141])

# calculate proportions of Tub+ and plot the results
fs = 20
counts2 = []
nobs2 = []
# colors to match Maheva's original plot
COL = np.array([[148, 225, 255, 255], [19, 129, 56, 255], [240, 224, 129, 255], [242, 242, 242, 255], [55, 36, 151, 255]])/255
plt.figure(figsize=(8, 10))
for i, loc in enumerate(locs):
    mea = np.mean(tub_bool[reg==loc])
    plt.bar(i, mea, ec='k', lw=2, fc=COL[i])
    plt.errorbar(i, mea, yerr=[[mea-conf_low[i]], [conf_high[i]-mea]], lw=2, color='k')
    count = int(np.sum(tub_bool[reg==loc]))
    counts2.append(count)
    nob = int(np.sum(reg==loc))
    nobs2.append(nob)
    plt.text(i, conf_high[i]+0.005, str(count)+'/'+str(nob), horizontalalignment='center', fontsize=fs)
    plt.text(i+0.05, mea + 0.002, 'ns', horizontalalignment='left', fontsize=fs)
#plt.xticks(range(5), locs)
plt.xticks([])
plt.yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25], ['0', '', '0.1', '', '0.2', ''], fontsize=fs)
# plt.yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35], ['0', '', '0.1', '', '0.2', '', '0.3', ''], fontsize=fs)
plt.ylabel('Proportion of Tubulin positive axons \n emerging through time', fontsize=fs)
plt.tight_layout()
plt.savefig(str(dirname)+'/tub_pos_axons.png')
plt.savefig('teste123')

print('one sided')
for i in range(5):
    print(np.sum(counts2[:i] + counts2[i+1 :])/np.sum(nobs2[:i] + nobs2[i+1 :]))
    print(proportions_ztest(counts2[i], nobs2[i], value=np.sum(counts2[:i] + counts2[i+1 :])/np.sum(nobs2[:i] + nobs[i+1 :])))

print('two sided')
for i in range(5):
    for j in range(5):
        if j > i:
            print(proportions_ztest([counts2[i], counts2[j]], [nobs2[i], nobs2[j]]))




