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
from is_notebook import is_notebook
from add_legend import add_legend

def assemblePerformanceTimeSeries():
    """Visualize performance time series for Case A of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script visualizes the performance time series quantities "
                    "as required by the CSP description."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-d', '--detailed', nargs='+', help='names of groups with detailed files', required=False)

    cmdArgs = vars(parser.parse_args())
    groups = set(cmdArgs["groups"])
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    detailed = []
    if cmdArgs["detailed"]:
        detailed = set(cmdArgs["detailed"])
        groups = groups.union(detailed)
    groups = sorted(list(groups))

    font = {'size' : 10, 'family': 'DejaVu Sans'}
    matplotlib.rc('font', **font)
    plt.rcParams['legend.title_fontsize'] = 'small'
    plt.rcParams['legend.fontsize'] = 'small'

    figT, axsT = plt.subplots(figsize=(5, 3))
    figF, axsF = plt.subplots(figsize=(5, 3))
    figM, axsM = plt.subplots(figsize=(5, 3))
    figD, axsD = plt.subplots(figsize=(5, 3))
    figN, axsN = plt.subplots(figsize=(5, 3))
    figR, axsR = plt.subplots(figsize=(5, 3))
    figL, axsL = plt.subplots(figsize=(5, 3))
    figRT, axsRT = plt.subplots(1, 2, figsize=(9, 3))
    figPub, axsPub = plt.subplots(2, 2, figsize=(9, 6))

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group.lower(), 'spe11a')
            if group.lower() in groups_and_colors:
                color = groups_and_colors[group.lower()]
            ls = '-'
            label = group
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1].lower(), 'spe11a', f'result{group[-1]}')
            if group[:-1].lower() in groups_and_colors:
                color = groups_and_colors[group[:-1].lower()]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'
            label = group[:-1]

        if group in detailed:
            fileName = os.path.join(baseFolder, 'spe11a_performance_time_series_detailed.csv')
        else:
            fileName = os.path.join(baseFolder, 'spe11a_performance_time_series.csv')
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60

        dtAvg = csvData[:, 1]
        dtAvg = np.convolve(dtAvg, [0.2, 0.2, 0.2, 0.2, 0.2], 'valid')
        dtAvg = np.insert(dtAvg, 0, csvData[0:2, 1])
        dtAvg = np.insert(dtAvg, -1, csvData[-2:, 1])
        axsT.plot(t, dtAvg, label=label, color=color, linestyle=ls)
        axsPub[0, 0].plot(t, dtAvg, label=label, color=color, linestyle=ls)
        axsF.plot(t, np.cumsum(csvData[:, 2]), label=label, color=color, linestyle=ls)
        axsM.plot(t, 1e3*csvData[:, 3], label=label, color=color, linestyle=ls)
        axsD.plot(t, csvData[:, 4], label=label, color=color, linestyle=ls)
        axsN.plot(t, np.cumsum(csvData[:, 5]), label=label, color=color, linestyle=ls)
        axsPub[1, 0].plot(t, np.cumsum(csvData[:, 5]), label=label, color=color, linestyle=ls)
        axsR.plot(t, np.cumsum(csvData[:, 6]), label=label, color=color, linestyle=ls)
        axsL.plot(t, np.cumsum(csvData[:, 7]), label=label, color=color, linestyle=ls)
        axsPub[1, 1].plot(t, np.cumsum(csvData[:, 7]), label=label, color=color, linestyle=ls)
        axsRT[0].plot(t, np.cumsum(csvData[:, 8]), label=label, color=color, linestyle=ls)
        axsPub[0, 1].plot(t, np.cumsum(csvData[:, 8]), label=label, color=color, linestyle=ls)
        axsRT[1].plot(t, np.cumsum(csvData[:, 9]), label=label, color=color, linestyle=ls)

    axsT.set_title(r'avg time step size')
    axsT.set_xlabel(r'time [h]')
    axsT.set_ylabel(r'step size [s]')
    axsT.set_xscale('log')
    axsT.set_yscale('log')
    axsT.set_xlim((9e-2, 2e2))
    axsT.set_ylim((1e-2, 1e3))
    add_legend(axsT)

    axsF.set_title(r'acc number of failed time steps')
    axsF.set_xlabel(r'time [h]')
    axsF.set_ylabel(r'failed steps [-]')
    axsF.set_xscale('log')
    axsF.set_yscale('log')
    axsF.set_xlim((9e-2, 2e2))
    add_legend(axsF)

    axsM.set_title(r'mass balance')
    axsM.set_xlabel(r'time [h]')
    axsM.set_ylabel(r'mass [g]')
    axsM.set_xscale('log')
    axsM.set_xlim((9e-2, 2e2))
    add_legend(axsM)

    axsD.set_title(r'avg degrees of freedom')
    axsD.set_xlabel(r'time [h]')
    axsD.set_ylabel(r'dof [-]')
    axsD.set_xscale('log')
    axsD.set_yscale('log')
    axsD.set_xlim((9e-2, 2e2))
    axsD.set_ylim((1e4, 1e7))
    add_legend(axsD)

    axsN.set_title(r'acc number of nonlinear iterations')
    axsN.set_xlabel(r'time [h]')
    axsN.set_ylabel(r'nonlinear iterations [-]')
    axsN.set_xscale('log')
    axsN.set_yscale('log')
    axsN.set_xlim((9e-2, 2e2))
    add_legend(axsN)

    axsR.set_title(r'acc number of local residual evaluations')
    axsR.set_xlabel(r'time [h]')
    axsR.set_ylabel(r'residual evaluations [-]')
    axsR.set_xscale('log')
    axsR.set_yscale('log')
    axsR.set_xlim((9e-2, 2e2))
    add_legend(axsR)

    axsL.set_title(r'acc number of linear iterations')
    axsL.set_xlabel(r'time [h]')
    axsL.set_ylabel(r'linear iterations [-]')
    axsL.set_xscale('log')
    axsL.set_yscale('log')
    axsL.set_xlim((9e-2, 2e2))
    add_legend(axsL)

    axsRT[0].set_title(r'acc runtime')
    axsRT[0].set_xlabel(r'time [h]')
    axsRT[0].set_ylabel(r'runtime [s]')
    axsRT[0].set_xscale('log')
    axsRT[0].set_yscale('log')
    axsRT[0].set_xlim((9e-2, 2e2))
    axsRT[1].set_title(r'acc time spent in linear solver')
    axsRT[1].set_xlabel(r'time [h]')
    axsRT[1].set_ylabel(r'runtime [s]')
    axsRT[1].set_xscale('log')
    axsRT[1].set_yscale('log')
    axsRT[1].set_xlim((9e-2, 2e2))
    add_legend(axsRT[1])

    axsPub[0, 0].set_title(r'avg time step size')
    axsPub[0, 0].set_ylabel(r'step size [s]')
    axsPub[0, 0].set_xscale('log')
    axsPub[0, 0].set_yscale('log')
    axsPub[0, 0].set_xlim([1e-1, 1.2e2])
    axsPub[0, 0].set_ylim([1e-2, 1e3])
    axsPub[0, 0].set_xticklabels([])
    axsPub[0, 1].set_title(r'acc runtime')
    axsPub[0, 1].set_ylabel(r'runtime [s]')
    axsPub[0, 1].set_xscale('log')
    axsPub[0, 1].set_yscale('log')
    axsPub[0, 1].set_xlim([1e-1, 1.2e2])
    axsPub[0, 1].yaxis.tick_right()
    axsPub[0, 1].yaxis.set_label_position('right')
    axsPub[0, 1].set_xticklabels([])
    axsPub[1, 0].set_title(r'acc number of nonlinear iterations')
    axsPub[1, 0].set_xlabel(r'time [h]')
    axsPub[1, 0].set_ylabel(r'nonlinear iterations [-]')
    axsPub[1, 0].set_xscale('log')
    axsPub[1, 0].set_yscale('log')
    axsPub[1, 0].set_xlim([1e-1, 1.2e2])
    axsPub[1, 1].set_title(r'acc number of linear iterations')
    axsPub[1, 1].set_xlabel(r'time [h]')
    axsPub[1, 1].set_ylabel(r'linear iterations [-]')
    axsPub[1, 1].set_xscale('log')
    axsPub[1, 1].set_yscale('log')
    axsPub[1, 1].set_xlim([1e-1, 1.2e2])
    axsPub[1, 1].yaxis.tick_right()
    axsPub[1, 1].yaxis.set_label_position('right')
    add_legend(axsPub[1, 1])

    figT.savefig(f'spe11a_time_series_tstep.png', bbox_inches='tight', dpi=300)
    figF.savefig(f'spe11a_time_series_fsteps.png', bbox_inches='tight', dpi=300)
    figM.savefig(f'spe11a_time_series_mass.png', bbox_inches='tight', dpi=300)
    figD.savefig(f'spe11a_time_series_dof.png', bbox_inches='tight', dpi=300)
    figN.savefig(f'spe11a_time_series_nliter.png', bbox_inches='tight', dpi=300)
    figR.savefig(f'spe11a_time_series_nres.png', bbox_inches='tight', dpi=300)
    figL.savefig(f'spe11a_time_series_liniter.png', bbox_inches='tight', dpi=300)
    figRT.savefig(f'spe11a_time_series_runtime.png', bbox_inches='tight', dpi=300)
    figPub.savefig(f'spe11a_performance_time_series.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    assemblePerformanceTimeSeries()
