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
    """Visualize time series for Case B of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script visualizes the time series quantities "
                    "as required by the CSP description."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups, taking reported numbers', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-c','--calculated', nargs='+', help='names of groups, taking calculated numbers for Boxes A and B')

    parser.add_argument('-t','--tablefolder', help='path to folder containing calculated tables')

    cmdArgs = vars(parser.parse_args())
    groups = [x.lower() for x in [x.lower() for x in cmdArgs["groups"]]]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    calculated = []
    if cmdArgs["calculated"]:
        calculated = [x.lower() for x in cmdArgs["calculated"]]
        groups = sorted(groups + calculated)
        tableFolder = cmdArgs["tablefolder"]

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    figP, axsP = plt.subplots(1, 2, figsize=(9, 3))
    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))
    figC, axsC = plt.subplots(figsize=(5, 3))
    figT, axsT = plt.subplots(1, 2, figsize=(9, 3))

    if calculated:
        mobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_mobA_from_spatial_maps.csv')
        mobileFromSpatialMapsA = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
        tSpatialMaps = mobileFromSpatialMapsA['time_s']/60/60/24/365
        immobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_immA_from_spatial_maps.csv')
        immobileFromSpatialMapsA = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
        dissolvedFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_dissA_from_spatial_maps.csv')
        dissolvedFromSpatialMapsA = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
        sealFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_sealA_from_spatial_maps.csv')
        sealFromSpatialMapsA = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

        mobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_mobB_from_spatial_maps.csv')
        mobileFromSpatialMapsB = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
        immobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_immB_from_spatial_maps.csv')
        immobileFromSpatialMapsB = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
        dissolvedFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_dissB_from_spatial_maps.csv')
        dissolvedFromSpatialMapsB = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
        sealFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_sealB_from_spatial_maps.csv')
        sealFromSpatialMapsB = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

        convectionFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_mC_from_spatial_maps.csv')
        convectionFromSpatialMaps = np.genfromtxt(convectionFromSpatialMapsFileName, delimiter=',', names=True)

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group, 'spe11b')
            if group in groups_and_colors:
                color = groups_and_colors[group]
            ls = '-'
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1], 'spe11b', f'result{group[-1]}')
            if group[:-1] in groups_and_colors:
                color = groups_and_colors[group[:-1]]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'

        fileName = os.path.join(baseFolder, 'spe11b_time_series.csv')
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365

        # scale pressure to bars
        if len(csvData[0]) > 1:
            axsP[0].plot(t, 1e-5*csvData[:, 1], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 2:
            axsP[1].plot(t, 1e-5*csvData[:, 2], label=group, color=color, linestyle=ls)

        # scale mass to kilotons
        if group in calculated:
            columnName = group.replace('-', '')
            axsA[0, 0].plot(tSpatialMaps, 1e-6*mobileFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsA[0, 1].plot(tSpatialMaps, 1e-6*immobileFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsA[1, 0].plot(tSpatialMaps, 1e-6*dissolvedFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsA[1, 1].plot(tSpatialMaps, 1e-6*sealFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
        else:
            axsA[0, 0].plot(t, 1e-6*csvData[:, 3], label=group, color=color, linestyle=ls)
            axsA[0, 1].plot(t, 1e-6*csvData[:, 4], label=group, color=color, linestyle=ls)
            axsA[1, 0].plot(t, 1e-6*csvData[:, 5], label=group, color=color, linestyle=ls)
            axsA[1, 1].plot(t, 1e-6*csvData[:, 6], label=group, color=color, linestyle=ls)
            # detect if immobile CO2 has been evaluated wrong potentially
            if max(1e-6*csvData[:, 4]) > 1:
                print(f"{group} potentially used inconsistent evaluation of immobile CO2.")

        if group in calculated:
            axsB[0, 0].plot(tSpatialMaps, 1e-6*mobileFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsB[0, 1].plot(tSpatialMaps, 1e-6*immobileFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsB[1, 0].plot(tSpatialMaps, 1e-6*dissolvedFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsB[1, 1].plot(tSpatialMaps, 1e-6*sealFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
        else:
            axsB[0, 0].plot(t, 1e-6*csvData[:, 7], label=group, color=color, linestyle=ls)
            axsB[0, 1].plot(t, 1e-6*csvData[:, 8], label=group, color=color, linestyle=ls)
            axsB[1, 0].plot(t, 1e-6*csvData[:, 9], label=group, color=color, linestyle=ls)
            axsB[1, 1].plot(t, 1e-6*csvData[:, 10], label=group, color=color, linestyle=ls)

        # scale length to kilometers
        if group in calculated:
            axsC.plot(tSpatialMaps, 1e-3*convectionFromSpatialMaps[columnName], label=group + r'$^*$', color=color, linestyle=ls)
        else:
            axsC.plot(t, 1e-3*csvData[:, 11], label=group, color=color, linestyle=ls)

        # scale mass to tons
        if len(csvData[0]) > 12:
            axsT[0].plot(t, 1e-3*csvData[:, 12], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 13:
            axsT[1].plot(t, 1e-3*csvData[:, 13], label=group, color=color, linestyle=ls)

    axsP[0].set_title(r'sensor 1')
    axsP[0].set_xlabel(r'time [y]')
    axsP[0].set_ylabel(r'pressure [bar]')
    axsP[0].set_xscale(r'log')
    axsP[0].set_xlim((1e0, 1e3))
    axsP[1].set_title(r'sensor 2')
    axsP[1].set_xlabel(r'time [y]')
    axsP[1].set_xscale(r'log')
    axsP[1].set_ylabel(r'pressure [bar]')
    axsP[1].yaxis.tick_right()
    axsP[1].yaxis.set_label_position('right')
    axsP[1].set_xlim((1e0, 1e3))
    handles, labels = axsP[1].get_legend_handles_labels()
    figP.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
    figP.tight_layout()
    figP.savefig('spe11b_time_series_pressure.png', bbox_inches='tight', dpi=120)

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO2')
    axsA[0, 0].set_ylabel(r'mass [kt]')
    axsA[0, 0].set_xticklabels([])
    axsA[0, 0].set_xscale(r'log')
    axsA[0, 0].set_xlim((1e0, 1e3))
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO2')
    axsA[0, 1].set_xticklabels([])
    axsA[0, 1].set_ylabel(r'mass [kt]')
    axsA[0, 1].yaxis.tick_right()
    axsA[0, 1].yaxis.set_label_position('right')
    axsA[0, 1].set_xscale(r'log')
    axsA[0, 1].set_xlim((1e0, 1e3))
    axsA[1, 0].set_title(r'Box A: dissolved CO2')
    axsA[1, 0].set_xlabel(r'time [y]')
    axsA[1, 0].set_ylabel(r'mass [kt]')
    axsA[1, 0].set_xscale(r'log')
    axsA[1, 0].set_xlim((1e0, 1e3))
    axsA[1, 1].set_title(r'Box A: CO2 in the seal facies')
    axsA[1, 1].set_xlabel(r'time [y]')
    axsA[1, 1].set_ylim(0, 5.0e-1)
    axsA[1, 1].set_ylabel(r'mass [kt]')
    axsA[1, 1].yaxis.tick_right()
    axsA[1, 1].yaxis.set_label_position('right')
    axsA[1, 1].set_xscale(r'log')
    axsA[1, 1].set_xlim((1e0, 1e3))
    if calculated:
        axsA[1, 1].plot([], [], ' ', label=r'$^*$from dense data')
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
    figA.tight_layout()
    figA.savefig('spe11b_time_series_boxA.png', bbox_inches='tight', dpi=120)

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
    axsB[0, 0].set_ylabel(r'mass [kt]')
    axsB[0, 0].set_xticklabels([])
    axsB[0, 0].set_xscale(r'log')
    axsB[0, 0].set_xlim((1e1, 1e3))
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
    axsB[0, 1].set_xticklabels([])
    axsB[0, 1].set_ylabel(r'mass [kt]')
    axsB[0, 1].yaxis.tick_right()
    axsB[0, 1].yaxis.set_label_position('right')
    axsB[0, 1].set_xscale(r'log')
    axsB[0, 1].set_xlim((1e1, 1e3))
    axsB[1, 0].set_title(r'Box B: dissolved CO2')
    axsB[1, 0].set_xlabel(r'time [y]')
    axsB[1, 0].set_ylabel(r'mass [kt]')
    axsB[1, 0].set_xscale(r'log')
    axsB[1, 0].set_xlim((1e1, 1e3))
    axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
    axsB[1, 1].set_xlabel(r'time [y]')
    axsB[1, 1].set_ylabel(r'mass [kt]')
    axsB[1, 1].yaxis.tick_right()
    axsB[1, 1].yaxis.set_label_position('right')
    axsB[1, 1].set_xscale(r'log')
    axsB[1, 1].set_xlim((1e1, 1e3))
    axsB[1, 1].set_ylim((0, 3.0e-1))
    if calculated:
        axsB[1, 1].plot([], [], ' ', label=r'$^*$from dense data')
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
    figB.tight_layout()
    figB.savefig('spe11b_time_series_boxB.png', bbox_inches='tight', dpi=120)

    axsC.set_title(r'Box C: convection')
    axsC.set_xlabel(r'time [y]')
    axsC.set_ylabel(r'$M$ [km]')
    axsC.set_xscale(r'log')
    axsC.set_xlim((1e1, 1e3))
    axsC.set_ylim((0, 20))
    if calculated:
        axsC.plot([], [], ' ', label=r'$^*$from dense data')
    axsC.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
    figC.savefig('spe11b_time_series_boxC.png', bbox_inches='tight', dpi=120)

    axsT[0].set_title(r'CO2 in sealing units')
    axsT[0].set_xlabel(r'time [y]')
    axsT[0].set_ylabel(r'mass [t]')
#    axsT[0].set_ylim(0, 1e3)
    axsT[0].set_xscale(r'log')
    axsT[0].set_xlim((1e0, 1e3))
    axsT[1].set_title(r'CO2 in boundary volumes')
    axsT[1].set_xlabel(r'time [y]')
    axsT[1].set_ylabel(r'mass [t]')
#    axsT[1].set_ylim(0, 1e2)
    axsT[1].set_xscale(r'log')
    axsT[1].yaxis.tick_right()
    axsT[1].yaxis.set_label_position('right')
    axsT[1].set_xlim((5e2, 1e3))
    handles, labels = axsT[1].get_legend_handles_labels()
    figT.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
    figT.tight_layout()
    figT.savefig('spe11b_time_series_seal.png', bbox_inches='tight', dpi=120)

if __name__ == "__main__":
    assembleTimeSeries()
