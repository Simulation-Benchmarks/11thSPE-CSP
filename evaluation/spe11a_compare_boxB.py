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

def compareConvection():
    """Compare different CO2 phase distributions in Box B for Case A of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script compares the reported values with the ones derived from spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-t','--tablefolder', help='path to folder containing calculated tables')

    cmdArgs = vars(parser.parse_args())
    groups = [x.lower() for x in cmdArgs["groups"]]
    folder = cmdArgs["folder"]
    tableFolder = cmdArgs["tablefolder"]

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    figM, axsM = plt.subplots(1, 2, figsize=(9, 3))
    figI, axsI = plt.subplots(1, 2, figsize=(9, 3))
    figD, axsD = plt.subplots(1, 2, figsize=(9, 3))
    figS, axsS = plt.subplots(1, 2, figsize=(9, 3))

    mobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11a_mobile_boxB_from_spatial_maps.csv')
    mobileFromSpatialMaps = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
    tSpatialMaps = mobileFromSpatialMaps['time_s']/60/60
    immobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11a_immobile_boxB_from_spatial_maps.csv')
    immobileFromSpatialMaps = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
    dissolvedFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11a_dissolved_boxB_from_spatial_maps.csv')
    dissolvedFromSpatialMaps = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
    sealFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11a_seal_boxB_from_spatial_maps.csv')
    sealFromSpatialMaps = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if not group[-1].isnumeric():
            baseFolder = os.path.join(folder, group, 'spe11a')
            if group in groups_and_colors:
                color = groups_and_colors[group]
            ls = '-'
        else:
            baseFolder = os.path.join(folder, group[:-1], 'spe11a', f'result{group[-1]}')
            if group[:-1] in groups_and_colors:
                color = groups_and_colors[group[:-1]]
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

        csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)
        t = csvData[:, 0]/60/60

        # scale mass to grams
        if len(csvData[0]) > 3:
            axsM[0].plot(t, 1e3*csvData[:, 7], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 4:
            axsI[0].plot(t, 1e3*csvData[:, 8], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 5:
            axsD[0].plot(t, 1e3*csvData[:, 9], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 6:
            axsS[0].plot(t, 1e3*csvData[:, 10], label=group, color=color, linestyle=ls)

        columnName = group.replace('-', '')
        axsM[1].plot(tSpatialMaps, 1e3*mobileFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)
        axsI[1].plot(tSpatialMaps, 1e3*immobileFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)
        axsD[1].plot(tSpatialMaps, 1e3*dissolvedFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)
        axsS[1].plot(tSpatialMaps, 1e3*sealFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)

    axsM[0].set_title(r'reported')
    axsM[0].set_ylabel(r'mass [g]')
    axsM[0].set_xscale(r'log')
    axsM[0].set_xlabel(r'time [h]')
    axsM[0].set_xlim((1e0, 1.2e2))
    axsM[0].set_ylim((-1e-2, 4.7e-1))
    axsM[1].set_title(r'from spatial maps')
    axsM[1].set_ylabel(r'mass [g]')
    axsM[1].set_xscale(r'log')
    axsM[1].set_xlabel(r'time [h]')
    axsM[1].set_xlim((1e0, 1.2e2))
    axsM[1].set_ylim((-1e-2, 4.7e-1))
    axsM[1].yaxis.tick_right()
    axsM[1].yaxis.set_label_position('right')
    figM.suptitle(r'Box B: mobile gaseous CO2')
    handles, labels = axsM[1].get_legend_handles_labels()
    figM.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figM.tight_layout()
    figM.savefig('spe11a_compare_mobile_boxB.png', bbox_inches='tight', dpi=300)

    axsI[0].set_title(r'reported')
    axsI[0].set_ylabel(r'mass [g]')
    axsI[0].set_xscale(r'log')
    axsI[0].set_xlabel(r'time [h]')
    axsI[0].set_xlim((1e0, 1.2e2))
    axsI[0].set_ylim((-1e-3, 4.0e-2))
    axsI[1].set_title(r'from spatial maps')
    axsI[1].set_ylabel(r'mass [g]')
    axsI[1].set_xscale(r'log')
    axsI[1].set_xlabel(r'time [h]')
    axsI[1].set_xlim((1e0, 1.2e2))
    axsI[1].set_ylim((-1e-3, 4.0e-2))
    axsI[1].yaxis.tick_right()
    axsI[1].yaxis.set_label_position('right')
    figI.suptitle(r'Box B: immobile gaseous CO2')
    handles, labels = axsI[1].get_legend_handles_labels()
    figI.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figI.tight_layout()
    figI.savefig('spe11a_compare_immobile_boxB.png', bbox_inches='tight', dpi=300)

    axsD[0].set_title(r'reported')
    axsD[0].set_ylabel(r'mass [g]')
    axsD[0].set_xscale(r'log')
    axsD[0].set_xlabel(r'time [h]')
    axsD[0].set_xlim((1e0, 1.2e2))
    axsD[0].set_ylim((-1e-2, 1.1e0))
    axsD[1].set_title(r'from spatial maps')
    axsD[1].set_ylabel(r'mass [g]')
    axsD[1].set_xscale(r'log')
    axsD[1].set_xlabel(r'time [h]')
    axsD[1].set_xlim((1e0, 1.2e2))
    axsD[1].set_ylim((-1e-2, 1.1e0))
    axsD[1].yaxis.tick_right()
    axsD[1].yaxis.set_label_position('right')
    figD.suptitle(r'Box B: dissolved CO2')
    handles, labels = axsD[1].get_legend_handles_labels()
    figD.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figD.tight_layout()
    figD.savefig('spe11a_compare_dissolved_boxB.png', bbox_inches='tight', dpi=300)

    axsS[0].set_title(r'reported')
    axsS[0].set_ylabel(r'mass [g]')
    axsS[0].set_xscale(r'log')
    axsS[0].set_xlabel(r'time [h]')
    axsS[0].set_xlim((1e0, 1.2e2))
    axsS[0].set_ylim((-1e-2, 2.5e-1))
    axsS[1].set_title(r'from spatial maps')
    axsS[1].set_ylabel(r'mass [g]')
    axsS[1].set_xscale(r'log')
    axsS[1].set_xlabel(r'time [h]')
    axsS[1].set_xlim((1e0, 1.2e2))
    axsS[1].set_ylim((-1e-2, 2.5e-1))
    axsS[1].yaxis.tick_right()
    axsS[1].yaxis.set_label_position('right')
    figS.suptitle(r'Box B: CO2 in the seal facies')
    handles, labels = axsS[1].get_legend_handles_labels()
    figS.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figS.tight_layout()
    figS.savefig('spe11a_compare_seal_boxB.png', bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    compareConvection()
