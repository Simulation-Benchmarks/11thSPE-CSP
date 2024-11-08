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
    """Compare different convection calculations for Case B of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script compares the different convection calculations."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    folder = cmdArgs["folder"]

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    fig, axs = plt.subplots(figsize=(5, 3))

    fromSpatialMapsFileName = os.path.join(folder, 'evaluation', 'spe11b_convection_from_spatial_maps.csv')
    fromSpatialMaps = np.genfromtxt(fromSpatialMapsFileName, delimiter=',', names=True)
    tSpatialMaps = fromSpatialMaps['time_s']/60/60/24/365

    typicalFileName = os.path.join(folder, 'ut-csee-pge', 'spe11b', 'result2', 'spe11b_time_series.csv')
    typical = np.genfromtxt(typicalFileName, delimiter=',')
    axs.plot(typical[:, 0]/60/60/24/365, 1e-3*typical[:, 11], label='typical', color='k', linestyle='-')

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if group[-2] != '-':
            baseFolder = os.path.join(folder, group.lower(), 'spe11b')
            if group.lower() in groups_and_colors:
                color = groups_and_colors[group.lower()]
        else:
            baseFolder = os.path.join(folder, group[:-2].lower(), 'spe11b', f'result{group[-1]}')
            if group[:-2].lower() in groups_and_colors:
                color = groups_and_colors[group[:-2].lower()]

        fileName = os.path.join(baseFolder, 'spe11b_time_series.csv')
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365

        # scale length to kilometers
        if len(csvData[0]) > 11:
            axs.plot(t, 1e-3*csvData[:, 11], label=group + ' reported', color=color, linestyle='-')

        columnName = group.lower().replace('-', '')
        axs.plot(tSpatialMaps, 1e-3*fromSpatialMaps[columnName], label=group + ' calculated', color=color, linestyle='--')

    axs.set_title(r'Box C: convection')
    axs.set_xlabel(r'time [y]')
    axs.set_ylabel(r'$M$ [km]')
    axs.set_xscale(r'log')
    axs.set_yscale(r'log')
    axs.set_xlim((1e0, 1e3))
    axs.set_ylim((1e-3, 1.2e1))
    axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig('spe11b_compare_convection.png', bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    compareConvection()
