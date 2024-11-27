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
    """Compare different CO2 phase distributions in Box A for Case B of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script compares the reported values with the ones derived from spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    folder = cmdArgs["folder"]

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    figM, axsM = plt.subplots(1, 2, figsize=(9, 3))
    figI, axsI = plt.subplots(1, 2, figsize=(9, 3))
    figD, axsD = plt.subplots(1, 2, figsize=(9, 3))
    figS, axsS = plt.subplots(1, 2, figsize=(9, 3))

    mobileFromSpatialMapsFileName = 'spe11b_mobile_boxA_from_spatial_maps.csv'
    mobileFromSpatialMaps = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
    tSpatialMaps = mobileFromSpatialMaps['time_s']/60/60/24/365
    immobileFromSpatialMapsFileName = 'spe11b_immobile_boxA_from_spatial_maps.csv'
    immobileFromSpatialMaps = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
    dissolvedFromSpatialMapsFileName = 'spe11b_dissolved_boxA_from_spatial_maps.csv'
    dissolvedFromSpatialMaps = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
    sealFromSpatialMapsFileName = 'spe11b_seal_boxA_from_spatial_maps.csv'
    sealFromSpatialMaps = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if group[-2] != '-':
            baseFolder = os.path.join(folder, group.lower(), 'spe11b')
            if group.lower() in groups_and_colors:
                color = groups_and_colors[group.lower()]
            ls = '-'
        else:
            baseFolder = os.path.join(folder, group[:-2].lower(), 'spe11b', f'result{group[-1]}')
            if group[:-2].lower() in groups_and_colors:
                color = groups_and_colors[group[:-2].lower()]
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

        csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365

        # scale mass to kilotons
        if len(csvData[0]) > 3:
            axsM[0].plot(t, 1e-6*csvData[:, 3], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 4:
            axsI[0].plot(t, 1e-6*csvData[:, 4], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 5:
            axsD[0].plot(t, 1e-6*csvData[:, 5], label=group, color=color, linestyle=ls)
        if len(csvData[0]) > 6:
            axsS[0].plot(t, 1e-6*csvData[:, 6], label=group, color=color, linestyle=ls)

        columnName = group.lower().replace('-', '')
        axsM[1].plot(tSpatialMaps, 1e-6*mobileFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)
        axsI[1].plot(tSpatialMaps, 1e-6*immobileFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)
        axsD[1].plot(tSpatialMaps, 1e-6*dissolvedFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)
        axsS[1].plot(tSpatialMaps, 1e-6*sealFromSpatialMaps[columnName], label=group, color=color, linestyle=ls)

    axsM[0].set_title(r'reported')
    axsM[0].set_ylabel(r'mass [kt]')
    axsM[0].set_xscale(r'log')
    axsM[0].set_xlim((1e0, 1e3))
    axsM[0].set_ylim((-1e0, 42))
    axsM[1].set_title(r'from spatial maps')
    axsM[1].set_ylabel(r'mass [kt]')
    axsM[1].set_xscale(r'log')
    axsM[1].set_xlim((1e0, 1e3))
    axsM[1].set_ylim((-1e0, 42))
    axsM[1].yaxis.tick_right()
    axsM[1].yaxis.set_label_position('right')
    figM.suptitle(r'Box A: mobile gaseous CO2')
    handles, labels = axsM[1].get_legend_handles_labels()
    figM.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figM.tight_layout()
    figM.savefig('spe11b_compare_mobile_boxA.png', bbox_inches='tight', dpi=300)

    axsI[0].set_title(r'reported')
    axsI[0].set_ylabel(r'mass [kt]')
    axsI[0].set_xscale(r'log')
    axsI[0].set_xlim((1e0, 1e3))
    axsI[0].set_ylim((-1e-2, 9e-1))
    axsI[1].set_title(r'from spatial maps')
    axsI[1].set_ylabel(r'mass [kt]')
    axsI[1].set_xscale(r'log')
    axsI[1].set_xlim((1e0, 1e3))
    axsI[1].set_ylim((-1e-2, 9e-1))
    axsI[1].yaxis.tick_right()
    axsI[1].yaxis.set_label_position('right')
    figI.suptitle(r'Box A: immobile gaseous CO2')
    handles, labels = axsI[1].get_legend_handles_labels()
    figI.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figI.tight_layout()
    figI.savefig('spe11b_compare_immobile_boxA.png', bbox_inches='tight', dpi=300)

    axsD[0].set_title(r'reported')
    axsD[0].set_ylabel(r'mass [kt]')
    axsD[0].set_xscale(r'log')
    axsD[0].set_xlim((1e0, 1e3))
    axsD[0].set_ylim((-5e-1, 1.7e1))
    axsD[1].set_title(r'from spatial maps')
    axsD[1].set_ylabel(r'mass [kt]')
    axsD[1].set_xscale(r'log')
    axsD[1].set_xlim((1e0, 1e3))
    axsD[1].set_ylim((-5e-1, 1.7e1))
    axsD[1].yaxis.tick_right()
    axsD[1].yaxis.set_label_position('right')
    figD.suptitle(r'Box A: dissolved CO2')
    handles, labels = axsD[1].get_legend_handles_labels()
    figD.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figD.tight_layout()
    figD.savefig('spe11b_compare_dissolved_boxA.png', bbox_inches='tight', dpi=300)

    axsS[0].set_title(r'reported')
    axsS[0].set_ylabel(r'mass [kt]')
    axsS[0].set_xscale(r'log')
    axsS[0].set_xlim((1e0, 1e3))
    axsS[0].set_ylim((-5e-2, 1.0e0))
    axsS[1].set_title(r'from spatial maps')
    axsS[1].set_ylabel(r'mass [kt]')
    axsS[1].set_xscale(r'log')
    axsS[1].set_xlim((1e0, 1e3))
    axsS[1].set_ylim((-5e-2, 1.0e0))
    axsS[1].yaxis.tick_right()
    axsS[1].yaxis.set_label_position('right')
    figS.suptitle(r'Box A: CO2 in the seal facies')
    handles, labels = axsS[1].get_legend_handles_labels()
    figS.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    figS.tight_layout()
    figS.savefig('spe11b_compare_seal_boxA.png', bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    compareConvection()