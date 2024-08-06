# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors

def assembleTimeSeries():
    """Visualize time series for Case A of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script visualizes the time series quantities "
                    "as required by the CSP description."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    figP, axsP = plt.subplots(1, 2, figsize=(9, 3))
    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))
    figC, axsC = plt.subplots(figsize=(5, 3))
    figT, axsT = plt.subplots(figsize=(5, 3))

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if groupFolders:
            baseFolder = groupFolders[i]

        if group[-2] != '-':
            if not groupFolders:
                baseFolder = os.path.join(folder, group.lower())
            if group.lower() in groups_and_colors:
                color = groups_and_colors[group.lower()]
            ls = '-'
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-2].lower(), f'result{group[-1]}')
            if group[:-2].lower() in groups_and_colors:
                color = groups_and_colors[group[:-2].lower()]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'

        fileName = os.path.join(baseFolder, 'spe11a_time_series.csv')
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60

        # scale pressure to bars
        axsP[0].plot(t, 1e-5*csvData[:, 1], label=group, color=color, linestyle=ls)
        axsP[1].plot(t, 1e-5*csvData[:, 2], label=group, color=color, linestyle=ls)

        # scale mass to grams
        axsA[0, 0].plot(t, 1e3*csvData[:, 3], label=group, color=color, linestyle=ls)
        axsA[0, 1].plot(t, 1e3*csvData[:, 4], label=group, color=color, linestyle=ls)
        axsA[1, 0].plot(t, 1e3*csvData[:, 5], label=group, color=color, linestyle=ls)
        axsA[1, 1].plot(t, 1e3*csvData[:, 6], label=group, color=color, linestyle=ls)
        # detect if immobile CO2 has been evaluated wrong potentially
        if max(1e3*csvData[:, 4]) > 0.05:
            print(f"{group} potentially used wrong evaluation of immobile CO2.")

        axsB[0, 0].plot(t, 1e3*csvData[:, 7], label=group, color=color, linestyle=ls)
        axsB[0, 1].plot(t, 1e3*csvData[:, 8], label=group, color=color, linestyle=ls)
        axsB[1, 0].plot(t, 1e3*csvData[:, 9], label=group, color=color, linestyle=ls)
        axsB[1, 1].plot(t, 1e3*csvData[:, 10], label=group, color=color, linestyle=ls)

        axsC.plot(t, csvData[:, 11], label=group, color=color, linestyle=ls)

        # scale mass to grams
        axsT.plot(t, 1e3*csvData[:, 12], label=group, color=color, linestyle=ls)

    axsP[0].set_title(r'sensor 1')
    axsP[0].set_xlabel(r'time [h]')
    axsP[0].set_ylabel(r'pressure [bar]')
    axsP[0].set_xlim(1e-1, 7260.0/60)
    axsP[0].set_ylim(1.13, 1.2)
    axsP[0].set_xscale(r'log')
    axsP[1].set_title(r'sensor 2')
    axsP[1].set_xlabel(r'time [h]')
    axsP[1].set_xlim(1e-1, 7260.0/60)
    axsP[1].set_ylim(1.09, 1.16)
    axsP[1].set_xscale(r'log')
    axsP[1].set_ylabel(r'pressure [bar]')
    axsP[1].yaxis.tick_right()
    axsP[1].yaxis.set_label_position('right')
    handles, labels = axsP[1].get_legend_handles_labels()
    figP.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figP.tight_layout()
    figP.savefig('spe11a_time_series_pressure.png', bbox_inches='tight', dpi=300)

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO2')
    axsA[0, 0].set_ylabel(r'mass [g]')
    axsA[0, 0].set_xlim(1e-1, 7260.0/60)
    axsA[0, 0].set_xticklabels([])
    axsA[0, 0].set_xscale(r'log')
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO2')
    axsA[0, 1].set_xlim(1e-1, 7260.0/60)
    axsA[0, 1].set_xticklabels([])
    axsA[0, 1].set_xscale(r'log')
    axsA[0, 1].set_ylabel(r'mass [g]')
    axsA[0, 1].yaxis.tick_right()
    axsA[0, 1].yaxis.set_label_position('right')
    axsA[1, 0].set_title(r'Box A: dissolved CO2')
    axsA[1, 0].set_xlabel(r'time [h]')
    axsA[1, 0].set_ylabel(r'mass [g]')
    axsA[1, 0].set_xlim(1e-1, 7260.0/60)
    axsA[1, 0].set_xscale(r'log')
    axsA[1, 1].set_xscale(r'log')
    axsA[1, 1].set_title(r'Box A: CO2 in the seal facies')
    axsA[1, 1].set_xlabel(r'time [h]')
    axsA[1, 1].set_xlim(1e-1, 7260.0/60)
    axsA[1, 1].set_ylabel(r'mass [g]')
    axsA[1, 1].yaxis.tick_right()
    axsA[1, 1].yaxis.set_label_position('right')
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figA.tight_layout()
    figA.savefig('spe11a_time_series_boxA.png', bbox_inches='tight', dpi=300)

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
    axsB[0, 0].set_ylabel(r'mass [g]')
    axsB[0, 0].set_xlim(2e0, 7260.0/60)
    axsB[0, 0].set_ylim(-0.01, 0.15)
    axsB[0, 0].set_xticklabels([])
    axsB[0, 0].set_xscale(r'log')
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
    axsB[0, 1].set_xlim(2e0, 7260.0/60)
    axsB[0, 1].set_xticklabels([])
    axsB[0, 1].set_xscale(r'log')
    axsB[0, 1].set_ylabel(r'mass [g]')
    axsB[0, 1].yaxis.tick_right()
    axsB[0, 1].yaxis.set_label_position('right')
    axsB[1, 0].set_title(r'Box B: dissolved CO2')
    axsB[1, 0].set_xlabel(r'time [h]')
    axsB[1, 0].set_ylabel(r'mass [g]')
    axsB[1, 0].set_xlim(2e0, 7260.0/60)
    axsB[1, 0].set_xscale(r'log')
    axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
    axsB[1, 1].set_xlabel(r'time [h]')
    axsB[1, 1].set_xlim(2e0, 7260.0/60)
    axsB[1, 1].set_xscale(r'log')
    axsB[1, 1].set_ylabel(r'mass [g]')
    axsB[1, 1].set_ylim(-0.01, 0.1)
    axsB[1, 1].yaxis.tick_right()
    axsB[1, 1].yaxis.set_label_position('right')
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figB.tight_layout()
    figB.savefig('spe11a_time_series_boxB.png', bbox_inches='tight', dpi=300)

    axsC.set_title(r'Box C: convection')
    axsC.set_xlabel(r'time [h]')
    axsC.set_ylabel(r'$M$ [m]')
    axsC.set_xlim(1e0, 7260.0/60)
    axsC.set_ylim(0, 3.5)
    axsC.set_xscale(r'log')
    axsC.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    figC.savefig('spe11a_time_series_boxC.png', bbox_inches='tight', dpi=300)

    axsT.set_title(r'CO2 in sealing units')
    axsT.set_xlabel(r'time [h]')
    axsT.set_ylabel(r'mass [g]')
    axsT.set_xlim(1e-1, 7260.0/60)
    axsT.set_ylim(0, 0.65)
    axsT.set_xscale(r'log')
    axsT.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    figT.savefig('spe11a_time_series_seal.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    assembleTimeSeries()
