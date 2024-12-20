# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors


font = {'size' : 12}
matplotlib.rc('font', **font)

parser = argparse.ArgumentParser(
    description="This calculates the convection integral ((17) in the description) "
                "based on given spatial maps."
)

parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

parser.add_argument('-t','--tablefolder', help='path to folder containing calculated tables')

cmdArgs = vars(parser.parse_args())
groups = [x.lower() for x in cmdArgs["groups"]]
tableFolder = cmdArgs["tablefolder"]

fromSpatialMapsFileName = os.path.join(tableFolder, 'spe11a_convection_from_spatial_maps.csv')
csvData = np.genfromtxt(fromSpatialMapsFileName, delimiter=',', skip_header=1)
t = csvData[:, 0]/60/60

fig, axs = plt.subplots(figsize=(5, 3))

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

    axs.plot(t, csvData[:, i+1], label=group, color=color, linestyle=ls)

axs.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
axs.set_title(r'Box C: convection from spatial maps')
axs.set_xlabel(r'time [h]')
axs.set_ylabel(r'$M$ [m]')
axs.set_xscale('log')
#axs.set_ylim((0, 14))
fig.savefig('spe11a_convection_from_spatial_maps.png', bbox_inches='tight', dpi=300)
