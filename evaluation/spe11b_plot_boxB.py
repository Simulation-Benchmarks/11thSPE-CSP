import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors


font = {'size' : 12}
matplotlib.rc('font', **font)

parser = argparse.ArgumentParser(
    description="This plots CO2 distribution in Box B for Case B based on calculated tables."
)

parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

cmdArgs = vars(parser.parse_args())
groups = cmdArgs["groups"]

mobileData = np.genfromtxt('spe11b_mobile_boxB_from_spatial_maps.csv', delimiter=',', skip_header=1)
immobileData = np.genfromtxt('spe11b_immobile_boxB_from_spatial_maps.csv', delimiter=',', skip_header=1)
dissolvedData = np.genfromtxt('spe11b_dissolved_boxB_from_spatial_maps.csv', delimiter=',', skip_header=1)
sealData = np.genfromtxt('spe11b_seal_boxB_from_spatial_maps.csv', delimiter=',', skip_header=1)
t = mobileData[:, 0]/60/60/24/365

figB, axsB = plt.subplots(2, 2, figsize=(9, 6))

for i, group in zip(range(len(groups)), groups):
    color = f'C{i}'

    if group[-2] != '-':
        if group.lower() in groups_and_colors:
            color = groups_and_colors[group.lower()]
        ls = '-'
    else:
        if group[:-2].lower() in groups_and_colors:
            color = groups_and_colors[group[:-2].lower()]
        if group[-1] == '1': ls = '-'
        elif group[-1] == '2': ls = '--'
        elif group[-1] == '3': ls = '-.'
        elif group[-1] == '4': ls = ':'

    # scale to kt
    axs[0,0].plot(t, 1e-6*mobileData[:, i+1], label=group, color=color, linestyle=ls)
    axs[0,1].plot(t, 1e-6*immobileData[:, i+1], label=group, color=color, linestyle=ls)
    axs[1,0].plot(t, 1e-6*dissolvedData[:, i+1], label=group, color=color, linestyle=ls)
    axs[1,1].plot(t, 1e-6*sealData[:, i+1], label=group, color=color, linestyle=ls)

axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
axsB[0, 0].set_ylabel(r'mass [kt]')
#axsB[0, 0].set_ylim(-0.5, 4e1)
axsB[0, 0].set_xticklabels([])
axsB[0, 0].set_xscale(r'log')
axsB[0, 0].set_xlim((1e0, 1e3))
axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
axsB[0, 1].set_xticklabels([])
axsB[0, 1].set_ylabel(r'mass [kt]')
axsB[0, 1].yaxis.tick_right()
axsB[0, 1].yaxis.set_label_position('right')
axsB[0, 1].set_xscale(r'log')
axsB[0, 1].set_xlim((1e0, 1e3))
axsB[1, 0].set_title(r'Box B: dissolved CO2')
axsB[1, 0].set_xlabel(r'time [y]')
axsB[1, 0].set_ylabel(r'mass [kt]')
#axsB[1, 0].set_ylim(0, 1e1)
axsB[1, 0].set_xscale(r'log')
axsB[1, 0].set_xlim((1e0, 1e3))
axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
axsB[1, 1].set_xlabel(r'time [y]')
#axsB[1, 1].set_ylim(0, 2e-1)
axsB[1, 1].set_ylabel(r'mass [kt]')
axsB[1, 1].yaxis.tick_right()
axsB[1, 1].yaxis.set_label_position('right')
axsB[1, 1].set_xscale(r'log')
axsB[1, 1].set_xlim((1e0, 1e3))
handles, labels = axsB[1][1].get_legend_handles_labels()
figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
figB.tight_layout()
figB.savefig('spe11b_time_series_boxB_from_spatial_maps.png', bbox_inches='tight', dpi=300)
