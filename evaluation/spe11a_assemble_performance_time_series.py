#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors

def assembleTimeSeries():
    """Visualize performance time series for Case A of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script visualizes the performance time series quantities "
                    "as required by the CSP description."
    )
    parser.add_argument('-d','--detailed', default='False', help='flag if detailed files should be considered')

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=True)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    folder = cmdArgs["folder"]
    detailed = cmdArgs["detailed"]

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

    for group in groups:
        if group[-2] != '-':
            if detailed == "False":
                fileName = os.path.join(folder, group.lower(), 'spe11a_performance_time_series.csv')
            else:
                fileName = os.path.join(folder, group.lower(), 'spe11a_performance_time_series_detailed.csv')
            color = groups_and_colors[group.lower()]
            ls = '-'
        else:
            if detailed == "False":
                fileName = os.path.join(folder, group[:-2].lower(), f'result{group[-1]}', 'spe11a_performance_time_series.csv')
            else:
                fileName = os.path.join(folder, group[:-2].lower(), f'result{group[-1]}', 'spe11a_performance_time_series_detailed.csv')
            color = groups_and_colors[group[:-2].lower()]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'

        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60

        axsT.plot(t, csvData[:, 1], label=group, color=color, linestyle=ls)
        axsF.plot(t, csvData[:, 2], label=group, color=color, linestyle=ls)
        axsM.plot(t, 1e3*csvData[:, 3], label=group, color=color, linestyle=ls)
        axsD.plot(t, csvData[:, 4], label=group, color=color, linestyle=ls)
        axsN.plot(t, csvData[:, 5], label=group, color=color, linestyle=ls)
        axsR.plot(t, csvData[:, 6], label=group, color=color, linestyle=ls)
        axsL.plot(t, csvData[:, 7], label=group, color=color, linestyle=ls)
        axsRT[0].plot(t, csvData[:, 8], label=group, color=color, linestyle=ls)
        axsRT[1].plot(t, csvData[:, 9], label=group, color=color, linestyle=ls)

    axsT.set_title(r'time step size')
    axsT.set_xlabel(r'time [h]')
    axsT.set_ylabel(r'step size [s]')
    axsT.set_xscale('log')
    axsT.set_yscale('log')
    axsT.set_xlim((9e-2, 2e2))
    axsT.set_ylim((1e-2, 1e3))
    axsT.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsF.set_title(r'number of failed time steps')
    axsF.set_xlabel(r'time [h]')
    axsF.set_ylabel(r'failed steps [-]')
    axsF.set_xscale('log')
    axsF.set_yscale('log')
    axsF.set_xlim((9e-2, 2e2))
    axsF.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsM.set_title(r'mass balance')
    axsM.set_xlabel(r'time [h]')
    axsM.set_ylabel(r'mass [g]')
    axsM.set_xscale('log')
    axsM.set_xlim((9e-2, 2e2))
    axsM.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsD.set_title(r'degrees of freedom')
    axsD.set_xlabel(r'time [h]')
    axsD.set_ylabel(r'dof [-]')
    axsD.set_xscale('log')
    axsD.set_yscale('log')
    axsD.set_xlim((9e-2, 2e2))
    axsD.set_ylim((1e4, 1e7))
    axsD.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsN.set_title(r'number of nonlinear iterations')
    axsN.set_xlabel(r'time [h]')
    axsN.set_ylabel(r'nonlinear iterations [-]')
    axsN.set_xscale('log')
    axsN.set_yscale('log')
    axsN.set_xlim((9e-2, 2e2))
    axsN.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsR.set_title(r'number of local residual evaluations')
    axsR.set_xlabel(r'time [h]')
    axsR.set_ylabel(r'residual evaluations [-]')
    axsR.set_xscale('log')
    axsR.set_yscale('log')
    axsR.set_xlim((9e-2, 2e2))
    axsR.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsL.set_title(r'number of linear iterations')
    axsL.set_xlabel(r'time [h]')
    axsL.set_ylabel(r'linear iterations [-]')
    axsL.set_xscale('log')
    axsL.set_yscale('log')
    axsL.set_xlim((9e-2, 2e2))
    axsL.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    axsRT[0].set_title(r'runtime')
    axsRT[0].set_xlabel(r'time [h]')
    axsRT[0].set_ylabel(r'runtime [s]')
    axsRT[0].set_xscale('log')
    axsRT[0].set_yscale('log')
    axsRT[0].set_xlim((9e-2, 2e2))
    axsRT[1].set_title(r'time spent in linear solver')
    axsRT[1].set_xlabel(r'time [h]')
    axsRT[1].set_ylabel(r'runtime [s]')
    axsRT[1].set_xscale('log')
    axsRT[1].set_yscale('log')
    axsRT[1].set_xlim((9e-2, 2e2))
    handles, labels = axsRT[1].get_legend_handles_labels()
    figRT.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figRT.tight_layout()

    if detailed == "False":
        figT.savefig('spe11a_time_series_tstep.png', bbox_inches='tight', dpi=300)
        figF.savefig('spe11a_time_series_fsteps.png', bbox_inches='tight', dpi=300)
        figM.savefig('spe11a_time_series_mass.png', bbox_inches='tight', dpi=300)
        figD.savefig('spe11a_time_series_dof.png', bbox_inches='tight', dpi=300)
        figN.savefig('spe11a_time_series_nliter.png', bbox_inches='tight', dpi=300)
        figR.savefig('spe11a_time_series_nres.png', bbox_inches='tight', dpi=300)
        figL.savefig('spe11a_time_series_liniter.png', bbox_inches='tight', dpi=300)
        figRT.savefig('spe11a_time_series_runtime.png', bbox_inches='tight', dpi=300)
    else:
        figT.savefig('spe11a_time_series_tstep_detailed.png', bbox_inches='tight', dpi=300)
        figF.savefig('spe11a_time_series_fsteps_detailed.png', bbox_inches='tight', dpi=300)
        figM.savefig('spe11a_time_series_mass_detailed.png', bbox_inches='tight', dpi=300)
        figD.savefig('spe11a_time_series_dof_detailed.png', bbox_inches='tight', dpi=300)
        figN.savefig('spe11a_time_series_nliter_detailed.png', bbox_inches='tight', dpi=300)
        figR.savefig('spe11a_time_series_nres_detailed.png', bbox_inches='tight', dpi=300)
        figL.savefig('spe11a_time_series_liniter_detailed.png', bbox_inches='tight', dpi=300)
        figRT.savefig('spe11a_time_series_runtime_detailed.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    assembleTimeSeries()
