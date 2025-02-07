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

def assemblePerformanceTimeSeries():
    """Visualize performance time series for Case C of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script visualizes the performance time series quantities "
                    "as required by the CSP description."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-d', '--detailed', required=False, help='set to true if detailed files should be considered', action=argparse.BooleanOptionalAction)
    parser.set_defaults(detailed=False)

    cmdArgs = vars(parser.parse_args())
    groups = [x.lower() for x in cmdArgs["groups"]]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    det = cmdArgs["detailed"]

    detailedString = ''
    if det:
        detailedString = '_detailed'
    csvName = f'spe11c_performance_time_series{detailedString}.csv'

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    figT, axsT = plt.subplots(figsize=(5, 3))
    figF, axsF = plt.subplots(figsize=(5, 3))
    figM, axsM = plt.subplots(figsize=(5, 3))
    figD, axsD = plt.subplots(figsize=(5, 3))
    figN, axsN = plt.subplots(figsize=(5, 3))
    figR, axsR = plt.subplots(figsize=(5, 3))
    figL, axsL = plt.subplots(figsize=(5, 3))
    figRT, axsRT = plt.subplots(1, 2, figsize=(9, 3))

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group, 'spe11c')
            if group in groups_and_colors:
                color = groups_and_colors[group]
            ls = '-'
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1], 'spe11c', f'result{group[-1]}')
            if group[:-1] in groups_and_colors:
                color = groups_and_colors[group[:-1]]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'

        fileName = os.path.join(baseFolder, csvName)
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365

        # scale time to days
        dtAvg = csvData[:, 1]
#        dtAvg = np.convolve(dtAvg, [0.2, 0.2, 0.2, 0.2, 0.2], 'valid')
#        dtAvg = np.insert(dtAvg, 0, csvData[0:2, 1])
#        dtAvg = np.insert(dtAvg, -1, csvData[-2:, 1])
        axsT.plot(t, dtAvg/60/60/24, label=group, color=color, linestyle=ls)
        axsF.plot(t, np.cumsum(csvData[:, 2]), label=group, color=color, linestyle=ls)
        # scale mass to kilotons
        axsM.plot(t, 1e-6*csvData[:, 3], label=group, color=color, linestyle=ls)
        axsD.plot(t, csvData[:, 4], label=group, color=color, linestyle=ls)
        axsN.plot(t, np.cumsum(csvData[:, 5]), label=group, color=color, linestyle=ls)
        axsR.plot(t, np.cumsum(csvData[:, 6]), label=group, color=color, linestyle=ls)
        axsL.plot(t, np.cumsum(csvData[:, 7]), label=group, color=color, linestyle=ls)
        axsRT[0].plot(t, np.cumsum(csvData[:, 8]), label=group, color=color, linestyle=ls)
        axsRT[1].plot(t, np.cumsum(csvData[:, 9]), label=group, color=color, linestyle=ls)

    axsT.set_title(r'avg time step size')
    axsT.set_xlabel(r'time [y]')
    axsT.set_ylabel(r'step size [d]')
    axsT.set_xscale('log')
    axsT.set_yscale('log')
    axsT.set_xlim([1e-1, 1e3])
    axsT.set_ylim([1e0, 2e3])
    axsT.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsF.set_title(r'acc number of failed time steps')
    axsF.set_xlabel(r'time [y]')
    axsF.set_ylabel(r'failed steps [-]')
    axsF.set_xscale('log')
    axsF.set_yscale('log')
    axsF.set_xlim([1e-1, 1e3])
    axsF.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsM.set_title(r'mass balance')
    axsM.set_xlabel(r'time [y]')
    axsM.set_ylabel(r'mass [kt]')
    axsM.set_xscale('log')
    axsM.set_xlim([1e-1, 1e3])
    axsM.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsD.set_title(r'avg degrees of freedom')
    axsD.set_xlabel(r'time [y]')
    axsD.set_ylabel(r'dof [-]')
    axsD.set_xscale('log')
    axsD.set_yscale('log')
    axsD.set_xlim([1e-1, 1e3])
    axsD.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsN.set_title(r'acc number of nonlinear iterations')
    axsN.set_xlabel(r'time [y]')
    axsN.set_ylabel(r'nonlinear iterations [-]')
    axsN.set_xscale('log')
    axsN.set_yscale('log')
    axsN.set_xlim([1e-1, 1e3])
    axsN.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsR.set_title(r'acc number of local residual evaluations')
    axsR.set_xlabel(r'time [y]')
    axsR.set_ylabel(r'residual evaluations [-]')
    axsR.set_xscale('log')
    axsR.set_yscale('log')
    axsR.set_xlim([1e-1, 1e3])
    axsR.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsL.set_title(r'acc number of linear iterations')
    axsL.set_xlabel(r'time [y]')
    axsL.set_ylabel(r'linear iterations [-]')
    axsL.set_xscale('log')
    axsL.set_yscale('log')
    axsL.set_xlim([1e-1, 1e3])
    axsL.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsRT[0].set_title(r'acc runtime')
    axsRT[0].set_xlabel(r'time [y]')
    axsRT[0].set_ylabel(r'runtime [s]')
    axsRT[0].set_xscale('log')
    axsRT[0].set_yscale('log')
    axsRT[0].set_xlim([1e-1, 1e3])
    axsRT[1].set_title(r'acc time spent in linear solver')
    axsRT[1].set_xlabel(r'time [y]')
    axsRT[1].set_ylabel(r'runtime [s]')
    axsRT[1].set_xscale('log')
    axsRT[1].set_yscale('log')
    axsRT[1].set_xlim([1e-1, 1e3])
    handles, labels = axsRT[1].get_legend_handles_labels()
    figRT.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figRT.tight_layout()

    figT.savefig(f'spe11c_time_series_tstep{detailedString}.png', bbox_inches='tight', dpi=300)
    figF.savefig(f'spe11c_time_series_fsteps{detailedString}.png', bbox_inches='tight', dpi=300)
    figM.savefig(f'spe11c_time_series_mass{detailedString}.png', bbox_inches='tight', dpi=300)
    figD.savefig(f'spe11c_time_series_dof{detailedString}.png', bbox_inches='tight', dpi=300)
    figN.savefig(f'spe11c_time_series_nliter{detailedString}.png', bbox_inches='tight', dpi=300)
    figR.savefig(f'spe11c_time_series_nres{detailedString}.png', bbox_inches='tight', dpi=300)
    figL.savefig(f'spe11c_time_series_liniter{detailedString}.png', bbox_inches='tight', dpi=300)
    figRT.savefig(f'spe11c_time_series_runtime{detailedString}.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    assemblePerformanceTimeSeries()
