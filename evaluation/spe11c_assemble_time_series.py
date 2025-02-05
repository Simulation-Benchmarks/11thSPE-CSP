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
    """Visualize time series for Case C of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script visualizes the time series quantities "
                    "as required by the CSP description."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-c','--calculated', nargs='+', help='names of groups, taking calculated numbers for Boxes A and B')

    parser.add_argument('-t','--tablefolder', help='path to folder containing calculated tables')

    parser.add_argument('-p','--percentiles', nargs='+', help='plot area between the two given percentiles')

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    calculated = []
    if cmdArgs["calculated"]:
        calculated = cmdArgs["calculated"]
        groups = sorted(groups + calculated)
        tableFolder = cmdArgs["tablefolder"]
    if cmdArgs["percentiles"]:
        lowerPercentile = int(cmdArgs["percentiles"][0])
        upperPercentile = int(cmdArgs["percentiles"][1])
        tableFolder = cmdArgs["tablefolder"]

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    figP, axsP = plt.subplots(1, 2, figsize=(9, 3))
    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))
    figC, axsC = plt.subplots(figsize=(5, 3))
    figT, axsT = plt.subplots(1, 2, figsize=(9, 3))

    if calculated:
        mobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_mobA_from_spatial_maps.csv')
        mobileFromSpatialMapsA = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
        tSpatialMaps = mobileFromSpatialMapsA['time_s']/60/60/24/365
        immobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_immA_from_spatial_maps.csv')
        immobileFromSpatialMapsA = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
        dissolvedFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_dissA_from_spatial_maps.csv')
        dissolvedFromSpatialMapsA = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
        sealFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_sealA_from_spatial_maps.csv')
        sealFromSpatialMapsA = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

        mobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_mobB_from_spatial_maps.csv')
        mobileFromSpatialMapsB = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
        immobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_immB_from_spatial_maps.csv')
        immobileFromSpatialMapsB = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
        dissolvedFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_dissB_from_spatial_maps.csv')
        dissolvedFromSpatialMapsB = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
        sealFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_sealB_from_spatial_maps.csv')
        sealFromSpatialMapsB = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

        convectionFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_mC_from_spatial_maps.csv')
        convectionFromSpatialMaps = np.genfromtxt(convectionFromSpatialMapsFileName, delimiter=',', names=True)

    if cmdArgs["percentiles"]:
        p1PercentilesFileName = os.path.join(tableFolder, 'spe11c_p1_percentiles.csv')
        p1Percentiles = np.genfromtxt(p1PercentilesFileName, delimiter=',', names=True)
        tPercentiles = p1Percentiles['time_s']/60/60/24/365
        axsP[0].fill_between(tPercentiles, 1e-5*p1Percentiles[f'P{lowerPercentile}_Pa'], 1e-5*p1Percentiles[f'P{upperPercentile}_Pa'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        p2PercentilesFileName = os.path.join(tableFolder, 'spe11c_p2_percentiles.csv')
        p2Percentiles = np.genfromtxt(p2PercentilesFileName, delimiter=',', names=True)
        axsP[1].fill_between(tPercentiles, 1e-5*p2Percentiles[f'P{lowerPercentile}_Pa'], 1e-5*p2Percentiles[f'P{upperPercentile}_Pa'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        mobAPercentilesFileName = os.path.join(tableFolder, 'spe11c_mobA_percentiles.csv')
        mobAPercentiles = np.genfromtxt(mobAPercentilesFileName, delimiter=',', names=True)
        axsA[0, 0].fill_between(tPercentiles, 1e-9*mobAPercentiles[f'P{lowerPercentile}_kg'], 1e-9*mobAPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        immAPercentilesFileName = os.path.join(tableFolder, 'spe11c_immA_percentiles.csv')
        immAPercentiles = np.genfromtxt(immAPercentilesFileName, delimiter=',', names=True)
        axsA[0, 1].fill_between(tPercentiles, 1e-9*immAPercentiles[f'P{lowerPercentile}_kg'], 1e-9*immAPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        dissAPercentilesFileName = os.path.join(tableFolder, 'spe11c_dissA_percentiles.csv')
        dissAPercentiles = np.genfromtxt(dissAPercentilesFileName, delimiter=',', names=True)
        axsA[1, 0].fill_between(tPercentiles, 1e-9*dissAPercentiles[f'P{lowerPercentile}_kg'], 1e-9*dissAPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        sealAPercentilesFileName = os.path.join(tableFolder, 'spe11c_sealA_percentiles.csv')
        sealAPercentiles = np.genfromtxt(sealAPercentilesFileName, delimiter=',', names=True)
        axsA[1, 1].fill_between(tPercentiles, 1e-9*sealAPercentiles[f'P{lowerPercentile}_kg'], 1e-9*sealAPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        mobBPercentilesFileName = os.path.join(tableFolder, 'spe11c_mobB_percentiles.csv')
        mobBPercentiles = np.genfromtxt(mobBPercentilesFileName, delimiter=',', names=True)
        axsB[0, 0].fill_between(tPercentiles, 1e-9*mobBPercentiles[f'P{lowerPercentile}_kg'], 1e-9*mobBPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        immBPercentilesFileName = os.path.join(tableFolder, 'spe11c_immB_percentiles.csv')
        immBPercentiles = np.genfromtxt(immBPercentilesFileName, delimiter=',', names=True)
        axsB[0, 1].fill_between(tPercentiles, 1e-9*immBPercentiles[f'P{lowerPercentile}_kg'], 1e-9*immBPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        dissBPercentilesFileName = os.path.join(tableFolder, 'spe11c_dissB_percentiles.csv')
        dissBPercentiles = np.genfromtxt(dissBPercentilesFileName, delimiter=',', names=True)
        axsB[1, 0].fill_between(tPercentiles, 1e-9*dissBPercentiles[f'P{lowerPercentile}_kg'], 1e-9*dissBPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        sealBPercentilesFileName = os.path.join(tableFolder, 'spe11c_sealB_percentiles.csv')
        sealBPercentiles = np.genfromtxt(sealBPercentilesFileName, delimiter=',', names=True)
        axsB[1, 1].fill_between(tPercentiles, 1e-9*sealBPercentiles[f'P{lowerPercentile}_kg'], 1e-9*sealBPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        mCPercentilesFileName = os.path.join(tableFolder, 'spe11c_mC_percentiles.csv')
        mCPercentiles = np.genfromtxt(mCPercentilesFileName, delimiter=',', names=True)
        axsC.fill_between(tPercentiles, 1e-6*mCPercentiles[f'P{lowerPercentile}_m'], 1e-6*mCPercentiles[f'P{upperPercentile}_m'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        sealTotPercentilesFileName = os.path.join(tableFolder, 'spe11c_sealTot_percentiles.csv')
        sealTotPercentiles = np.genfromtxt(sealTotPercentilesFileName, delimiter=',', names=True)
        axsT[0].fill_between(tPercentiles, 1e-9*sealTotPercentiles[f'P{lowerPercentile}_kg'], 1e-9*sealTotPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')
        boundaryPercentilesFileName = os.path.join(tableFolder, 'spe11c_boundary_percentiles.csv')
        boundaryPercentiles = np.genfromtxt(boundaryPercentilesFileName, delimiter=',', names=True)
        axsT[1].fill_between(tPercentiles, 1e-9*boundaryPercentiles[f'P{lowerPercentile}_kg'], 1e-9*boundaryPercentiles[f'P{upperPercentile}_kg'], color='lightgray', label=f'[P{lowerPercentile}, P{upperPercentile}]')

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group.lower(), 'spe11c')
            if group.lower() in groups_and_colors:
                color = groups_and_colors[group.lower()]
            ls = '-'
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1].lower(), 'spe11c', f'result{group[-1]}')
            if group[:-1].lower() in groups_and_colors:
                color = groups_and_colors[group[:-1].lower()]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'

        fileName = os.path.join(baseFolder, 'spe11c_time_series.csv')
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365

        # scale pressure to bars
        if np.isnan(csvData[:, 1]).all():
            print(f'{group} only reported nan for p1.')
        elif all(val < 0 for val in csvData[:, 1]):
            print(f'{group} only reported negative values for p1.')
        else:
            axsP[0].plot(t, 1e-5*csvData[:, 1], label=group, color=color, linestyle=ls)
        if np.isnan(csvData[:, 2]).all():
            print(f'{group} only reported nan for p2.')
        elif all(val < 0 for val in csvData[:, 2]):
            print(f'{group} only reported negative values for p2.')
        else:
            axsP[1].plot(t, 1e-5*csvData[:, 2], label=group, color=color, linestyle=ls)

        # scale mass to megatons
        if group in calculated:
            columnName = group.lower().replace('-', '')
            axsA[0, 0].plot(tSpatialMaps, 1e-9*mobileFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsA[0, 1].plot(tSpatialMaps, 1e-9*immobileFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsA[1, 0].plot(tSpatialMaps, 1e-9*dissolvedFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsA[1, 1].plot(tSpatialMaps, 1e-9*sealFromSpatialMapsA[columnName], label=group + r'$^*$', color=color, linestyle=ls)
        else:
            if np.isnan(csvData[:, 3]).all():
                print(f'{group} only reported nan for mobA.')
            elif all(val < 0 for val in csvData[:, 3]):
                print(f'{group} only reported negative values for mobA.')
            else:
                axsA[0, 0].plot(t, 1e-9*csvData[:, 3], label=group, color=color, linestyle=ls)
            if np.isnan(csvData[:, 4]).all():
                print(f'{group} only reported nan for immA.')
            elif all(val < 0 for val in csvData[:, 4]):
                print(f'{group} only reported negative values for immA.')
            else:
                axsA[0, 1].plot(t, 1e-9*csvData[:, 4], label=group, color=color, linestyle=ls)
            if np.isnan(csvData[:, 5]).all():
                print(f'{group} only reported nan for dissA.')
            elif all(val < 0 for val in csvData[:, 5]):
                print(f'{group} only reported negative values for dissA.')
            else:
                axsA[1, 0].plot(t, 1e-9*csvData[:, 5], label=group, color=color, linestyle=ls)
            if np.isnan(csvData[:, 6]).all():
                print(f'{group} only reported nan for sealA.')
            elif all(val < 0 for val in csvData[:, 6]):
                print(f'{group} only reported negative values for sealA.')
            else:
                axsA[1, 1].plot(t, 1e-9*csvData[:, 6], label=group, color=color, linestyle=ls)
            # detect if immobile CO2 has been evaluated wrong potentially
            if max(1e-9*csvData[:, 4]) > 2:
                print(f"{group} potentially used inconsistent evaluation of immobile CO2.")

        if group in calculated:
            axsB[0, 0].plot(tSpatialMaps, 1e-9*mobileFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsB[0, 1].plot(tSpatialMaps, 1e-9*immobileFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsB[1, 0].plot(tSpatialMaps, 1e-9*dissolvedFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
            axsB[1, 1].plot(tSpatialMaps, 1e-9*sealFromSpatialMapsB[columnName], label=group + r'$^*$', color=color, linestyle=ls)
        else:
            if np.isnan(csvData[:, 7]).all():
                print(f'{group} only reported nan for mobB.')
            elif all(val < 0 for val in csvData[:, 7]):
                print(f'{group} only reported negative values for mobB.')
            else:
                axsB[0, 0].plot(t, 1e-9*csvData[:, 7], label=group, color=color, linestyle=ls)
            if np.isnan(csvData[:, 8]).all():
                print(f'{group} only reported nan for immB.')
            elif all(val < 0 for val in csvData[:, 8]):
                print(f'{group} only reported negative values for immB.')
            else:
                axsB[0, 1].plot(t, 1e-9*csvData[:, 8], label=group, color=color, linestyle=ls)
            if np.isnan(csvData[:, 9]).all():
                print(f'{group} only reported nan for dissB.')
            elif all(val < 0 for val in csvData[:, 9]):
                print(f'{group} only reported negative values for dissB.')
            else:
                axsB[1, 0].plot(t, 1e-9*csvData[:, 9], label=group, color=color, linestyle=ls)
            if np.isnan(csvData[:, 10]).all():
                print(f'{group} only reported nan for sealB.')
            elif all(val < 0 for val in csvData[:, 10]):
                print(f'{group} only reported negative values for sealB.')
            else:
                axsB[1, 1].plot(t, 1e-9*csvData[:, 10], label=group, color=color, linestyle=ls)

        # scale area to square km
        if group in calculated:
            axsC.plot(tSpatialMaps, 1e-6*convectionFromSpatialMaps[columnName], label=group + r'$^*$', color=color, linestyle=ls)
        else:
            if np.isnan(csvData[:, 11]).all():
                print(f'{group} only reported nan for mC.')
            elif all(val < 0 for val in csvData[:, 11]):
                print(f'{group} only reported negative values for mC.')
            else:
                axsC.plot(t, 1e-6*csvData[:, 11], label=group, color=color, linestyle=ls)

        # scale mass to megatons
        if np.isnan(csvData[:, 12]).all():
            print(f'{group} only reported nan for sealTot.')
        elif all(val < 0 for val in csvData[:, 12]):
            print(f'{group} only reported negative values for sealTot.')
        else:
            axsT[0].plot(t, 1e-9*csvData[:, 12], label=group, color=color, linestyle=ls)
        if np.isnan(csvData[:, 13]).all():
            print(f'{group} only reported nan for boundaryCO2.')
        elif all(val < 0 for val in csvData[:, 13]):
            print(f'{group} only reported negative values for boundaryCO2.')
        else:
            axsT[1].plot(t, 1e-9*csvData[:, 13], label=group, color=color, linestyle=ls)

    axsP[0].set_title(r'sensor 1')
    axsP[0].set_xlabel(r'time [y]')
    axsP[0].set_ylabel(r'pressure [bar]')
    axsP[0].set_xscale('log')
    axsP[0].set_xlim((1e0, 1e3))
    axsP[0].set_ylim((200, 320))
    axsP[1].set_title(r'sensor 2')
    axsP[1].set_xlabel(r'time [y]')
    axsP[1].set_xscale('log')
    axsP[1].set_xlim((1e0, 1e3))
    axsP[1].set_ylim((200, 280))
    axsP[1].set_ylabel(r'pressure [bar]')
    axsP[1].yaxis.tick_right()
    axsP[1].yaxis.set_label_position('right')
    handles, labels = axsP[1].get_legend_handles_labels()
    figP.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figP.tight_layout()
    figP.savefig('spe11c_time_series_pressure.png', bbox_inches='tight', dpi=300)

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO2')
    axsA[0, 0].set_ylabel(r'mass [Mt]')
    axsA[0, 0].set_xticklabels([])
    axsA[0, 0].set_xscale('log')
    axsA[0, 0].set_xlim((8e0, 1e3))
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO2')
    axsA[0, 1].set_xticklabels([])
    axsA[0, 1].set_xscale('log')
    axsA[0, 1].set_xlim((8e0, 1e3))
    axsA[0, 1].set_ylabel(r'mass [Mt]')
    axsA[0, 1].yaxis.tick_right()
    axsA[0, 1].yaxis.set_label_position('right')
    axsA[1, 0].set_title(r'Box A: dissolved CO2')
    axsA[1, 0].set_xlabel(r'time [y]')
    axsA[1, 0].set_ylabel(r'mass [Mt]')
    axsA[1, 0].set_xscale('log')
    axsA[1, 0].set_xlim((8e0, 1e3))
    axsA[1, 1].set_title(r'Box A: CO2 in the seal facies')
    axsA[1, 1].set_xlabel(r'time [y]')
    axsA[1, 1].set_xscale('log')
    axsA[1, 1].set_xlim((8e0, 1e3))
    axsA[1, 1].set_ylabel(r'mass [Mt]')
    axsA[1, 1].yaxis.tick_right()
    axsA[1, 1].yaxis.set_label_position('right')
    if calculated:
        axsA[1, 1].plot([], [], ' ', label=r'$^*$from dense data')
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figA.tight_layout()
    figA.savefig('spe11c_time_series_boxA.png', bbox_inches='tight', dpi=300)

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
    axsB[0, 0].set_ylabel(r'mass [Mt]')
    axsB[0, 0].set_xticklabels([])
    axsB[0, 0].set_xscale('log')
    axsB[0, 0].set_xlim((1e1, 1e3))
    axsB[0, 0].set_ylim((-0.002, 0.05))
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
    axsB[0, 1].set_xticklabels([])
    axsB[0, 1].set_xscale('log')
    axsB[0, 1].set_xlim((1e1, 1e3))
    axsB[0, 1].set_ylim((-0.002, 0.05))
    axsB[0, 1].set_ylabel(r'mass [Mt]')
    axsB[0, 1].yaxis.tick_right()
    axsB[0, 1].yaxis.set_label_position('right')
    axsB[1, 0].set_title(r'Box B: dissolved CO2')
    axsB[1, 0].set_xlabel(r'time [y]')
    axsB[1, 0].set_ylabel(r'mass [Mt]')
    axsB[1, 0].set_xscale('log')
    axsB[1, 0].set_xlim((1e1, 1e3))
    axsB[1, 0].set_ylim((-0.002, 0.05))
    axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
    axsB[1, 1].set_xlabel(r'time [y]')
    axsB[1, 1].set_xscale('log')
    axsB[1, 1].set_xlim((1e1, 1e3))
    axsB[1, 1].set_ylim((-1e-5, 1e-4))
    axsB[1, 1].set_ylabel(r'mass [Mt]')
    axsB[1, 1].yaxis.tick_right()
    axsB[1, 1].yaxis.set_label_position('right')
    if calculated:
        axsB[1, 1].plot([], [], ' ', label=r'$^*$from dense data')
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figB.tight_layout()
    figB.savefig('spe11c_time_series_boxB.png', bbox_inches='tight', dpi=300)

    axsC.set_title(r'Box C: convection')
    axsC.set_xlabel(r'time [y]')
    axsC.set_ylabel(r'$M$ [km$^2$]')
    axsC.set_xscale('log')
    axsC.set_xlim((8e1, 1e3))
    axsC.set_ylim((-1, 20))
    if calculated:
        axsC.plot([], [], ' ', label=r'$^*$from dense data')
    axsC.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    figC.savefig('spe11c_time_series_boxC.png', bbox_inches='tight', dpi=300)

    axsT[0].set_title(r'CO2 in sealing units')
    axsT[0].set_xlabel(r'time [y]')
    axsT[0].set_ylabel(r'mass [Mt]')
    axsT[0].set_xscale('log')
    axsT[0].set_xlim((1e0, 1e3))
    axsT[0].set_ylim((-0.1, 2))
    axsT[1].set_title(r'CO2 in boundary volumes')
    axsT[1].set_xlabel(r'time [y]')
    axsT[1].set_ylabel(r'mass [Mt]')
    axsT[1].set_xscale('log')
    axsT[1].set_xlim((1e2, 1e3))
    axsT[1].yaxis.tick_right()
    axsT[1].yaxis.set_label_position('right')
    handles, labels = axsT[1].get_legend_handles_labels()
    figT.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figT.tight_layout()
    figT.savefig('spe11c_time_series_seal.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    assembleTimeSeries()
