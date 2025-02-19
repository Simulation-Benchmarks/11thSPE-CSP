# SPDX-FileCopyrightText: 2025 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors


font = {'size' : 12}
matplotlib.rc('font', **font)

parser = argparse.ArgumentParser(
    description="This plots CO2 distribution in Box A for Case B based on calculated tables."
)

parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

cmdArgs = vars(parser.parse_args())
groups = [x.lower() for x in cmdArgs["groups"]]

mobileData = np.genfromtxt('spe11b_mobA_from_spatial_maps.csv', delimiter=',', skip_header=1)
immobileData = np.genfromtxt('spe11b_immA_from_spatial_maps.csv', delimiter=',', skip_header=1)
dissolvedData = np.genfromtxt('spe11b_dissA_from_spatial_maps.csv', delimiter=',', skip_header=1)
sealData = np.genfromtxt('spe11b_sealA_from_spatial_maps.csv', delimiter=',', skip_header=1)
t = mobileData[:, 0]/60/60/24/365

figA, axsA = plt.subplots(2, 2, figsize=(9, 6))

for i, group in zip(range(len(groups)), groups):
    color = f'C{i}'

    if not group[-1].isnumeric():
        if group in groups_and_colors:
            color = groups_and_colors[group]
        ls = '-'
    else:
        if group[:-1] in groups_and_colors:
            color = groups_and_colors[group[:-1]]
        if group[-1] == '1': ls = '-'
        elif group[-1] == '2': ls = '--'
        elif group[-1] == '3': ls = '-.'
        elif group[-1] == '4': ls = ':'

    # scale to kt
    axs[0,0].plot(t, 1e-6*mobileData[:, i+1], label=group, color=color, linestyle=ls)
    axs[0,1].plot(t, 1e-6*immobileData[:, i+1], label=group, color=color, linestyle=ls)
    axs[1,0].plot(t, 1e-6*dissolvedData[:, i+1], label=group, color=color, linestyle=ls)
    axs[1,1].plot(t, 1e-6*sealData[:, i+1], label=group, color=color, linestyle=ls)

axsA[0, 0].set_title(r'Box A: mobile gaseous CO$_2$')
axsA[0, 0].set_ylabel(r'mass [kt]')
#axsA[0, 0].set_ylim(-0.5, 4e1)
axsA[0, 0].set_xticklabels([])
axsA[0, 0].set_xscale(r'log')
axsA[0, 0].set_xlim((1e0, 1e3))
axsA[0, 1].set_title(r'Box A: immobile gaseous CO$_2$')
axsA[0, 1].set_xticklabels([])
axsA[0, 1].set_ylabel(r'mass [kt]')
axsA[0, 1].yaxis.tick_right()
axsA[0, 1].yaxis.set_label_position('right')
axsA[0, 1].set_xscale(r'log')
axsA[0, 1].set_xlim((1e0, 1e3))
axsA[1, 0].set_title(r'Box A: dissolved CO$_2$')
axsA[1, 0].set_xlabel(r'time [y]')
axsA[1, 0].set_ylabel(r'mass [kt]')
#axsA[1, 0].set_ylim(0, 1e1)
axsA[1, 0].set_xscale(r'log')
axsA[1, 0].set_xlim((1e0, 1e3))
axsA[1, 1].set_title(r'Box A: CO$_2$ in the seal facies')
axsA[1, 1].set_xlabel(r'time [y]')
#axsA[1, 1].set_ylim(0, 2e-1)
axsA[1, 1].set_ylabel(r'mass [kt]')
axsA[1, 1].yaxis.tick_right()
axsA[1, 1].yaxis.set_label_position('right')
axsA[1, 1].set_xscale(r'log')
axsA[1, 1].set_xlim((1e0, 1e3))
handles, labels = axsA[1][1].get_legend_handles_labels()
figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
figA.tight_layout()
figA.savefig('spe11b_time_series_boxA_from_spatial_maps.png', bbox_inches='tight', dpi=300)
