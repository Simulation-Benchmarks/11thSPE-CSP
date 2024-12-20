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
    """Compare different convection calculations for Case C of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script compares the different convection calculations."
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

    fig, axs = plt.subplots(figsize=(5, 3))

    fromSpatialMapsFileName = os.path.join(tableFolder, 'spe11c_convection_from_spatial_maps.csv')
    fromSpatialMaps = np.genfromtxt(fromSpatialMapsFileName, delimiter=',', names=True)
    tSpatialMaps = fromSpatialMaps['time_s']/60/60/24/365

    typicalFileName = os.path.join(folder, 'ut-csee-pge', 'spe11c', 'spe11c_time_series.csv')
    typical = np.genfromtxt(typicalFileName, delimiter=',')
    axs.plot(typical[:, 0]/60/60/24/365, 1e-6*typical[:, 11], label='typical', color='k', linestyle='-')

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if not group[-1].isnumeric():
            baseFolder = os.path.join(folder, group, 'spe11c')
            if group in groups_and_colors:
                color = groups_and_colors[group]
        else:
            baseFolder = os.path.join(folder, group[:-1], 'spe11c', f'result{group[-1]}')
            if group[:-1] in groups_and_colors:
                color = groups_and_colors[group[:-1]]

        fileName = os.path.join(baseFolder, 'spe11c_time_series.csv')
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365

        # scale length to square kilometers
        if len(csvData[0]) > 11:
            axs.plot(t, 1e-6*csvData[:, 11], label=group + ' reported', color=color, linestyle='-')

        columnName = group.replace('-', '')
        axs.plot(tSpatialMaps, 1e-6*fromSpatialMaps[columnName], label=group + ' calculated', color=color, linestyle='--')

    axs.set_title(r'Box C: convection')
    axs.set_xlabel(r'time [y]')
    axs.set_ylabel(r'$M$ [km$^2$]')
    axs.set_xscale(r'log')
    axs.set_xlim((5e0, 1e3))
    axs.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncols=2)
    fig.savefig('spe11c_compare_convection.png', bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    compareConvection()
