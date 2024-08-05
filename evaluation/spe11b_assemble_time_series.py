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

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=True)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    folder = cmdArgs["folder"]

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    figP, axsP = plt.subplots(1, 2, figsize=(9, 3))
    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))
    figC, axsC = plt.subplots(figsize=(5, 3))
    figT, axsT = plt.subplots(1, 2, figsize=(9, 3))

    for group in groups:
        if group[-2] != '-':
            fileName = folder + "/" + group.lower() + "/spe11b_time_series.csv"
            color = groups_and_colors[group.lower()]
            ls = '-'
        else:
            fileName = folder + "/" + group[:-2].lower() + "/result" + group[-1] + "/spe11b_time_series.csv"
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
        t = csvData[:, 0]/60/60/24/365

        # scale pressure to bars
        axsP[0].plot(t, 1e-5*csvData[:, 1], label=group, color=color, linestyle=ls)
        axsP[1].plot(t, 1e-5*csvData[:, 2], label=group, color=color, linestyle=ls)

        # scale mass to kilotons
        axsA[0, 0].plot(t, 1e-6*csvData[:, 3], label=group, color=color, linestyle=ls)
        axsA[0, 1].plot(t, 1e-6*csvData[:, 4], label=group, color=color, linestyle=ls)
        axsA[1, 0].plot(t, 1e-6*csvData[:, 5], label=group, color=color, linestyle=ls)
        axsA[1, 1].plot(t, 1e-6*csvData[:, 6], label=group, color=color, linestyle=ls)

        axsB[0, 0].plot(t, 1e-6*csvData[:, 7], label=group, color=color, linestyle=ls)
        axsB[0, 1].plot(t, 1e-6*csvData[:, 8], label=group, color=color, linestyle=ls)
        axsB[1, 0].plot(t, 1e-6*csvData[:, 9], label=group, color=color, linestyle=ls)
        axsB[1, 1].plot(t, 1e-6*csvData[:, 10], label=group, color=color, linestyle=ls)

        # scale length to meters
        axsC.plot(t, 1e-3*csvData[:, 11], label=group, color=color, linestyle=ls)

        # scale mass to tons
        axsT[0].plot(t, 1e-3*csvData[:, 12], label=group, color=color, linestyle=ls)
        axsT[1].plot(t, 1e-3*csvData[:, 13], label=group, color=color, linestyle=ls)

    axsP[0].set_title(r'sensor 1')
    axsP[0].set_xlabel(r'time [y]')
    axsP[0].set_ylabel(r'pressure [bar]')
    axsP[0].set_xscale(r'log')
    axsP[0].set_xlim((1e0, 1e3))
    axsP[0].set_ylim((270, 450))
    axsP[1].set_title(r'sensor 2')
    axsP[1].set_xlabel(r'time [y]')
    axsP[1].set_xscale(r'log')
    axsP[1].set_ylabel(r'pressure [bar]')
    axsP[1].yaxis.tick_right()
    axsP[1].yaxis.set_label_position('right')
    axsP[1].set_xlim((1e0, 1e3))
    axsP[1].set_ylim((210, 400))
    handles, labels = axsP[1].get_legend_handles_labels()
    figP.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figP.tight_layout()
    figP.savefig('spe11b_time_series_pressure.png', bbox_inches='tight', dpi=300)

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO2')
    axsA[0, 0].set_ylabel(r'mass [kt]')
    axsA[0, 0].set_ylim(0, 5e1)
    axsA[0, 0].set_xticklabels([])
    axsA[0, 0].set_xscale(r'log')
    axsA[0, 0].set_xlim((1e0, 1e3))
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO2')
    axsA[0, 1].set_ylim(0, 6e0)
    axsA[0, 1].set_xticklabels([])
    axsA[0, 1].set_ylabel(r'mass [kt]')
    axsA[0, 1].yaxis.tick_right()
    axsA[0, 1].yaxis.set_label_position('right')
    axsA[0, 1].set_xscale(r'log')
    axsA[0, 1].set_xlim((1e0, 1e3))
    axsA[1, 0].set_title(r'Box A: dissolved CO2')
    axsA[1, 0].set_xlabel(r'time [y]')
    axsA[1, 0].set_ylabel(r'mass [kt]')
    axsA[1, 0].set_ylim(0, 1e1)
    axsA[1, 0].set_xscale(r'log')
    axsA[1, 0].set_xlim((1e0, 1e3))
    axsA[1, 1].set_title(r'Box A: CO2 in the seal facies')
    axsA[1, 1].set_xlabel(r'time [y]')
    axsA[1, 1].set_ylim(0, 2e-1)
    axsA[1, 1].set_ylabel(r'mass [kt]')
    axsA[1, 1].yaxis.tick_right()
    axsA[1, 1].yaxis.set_label_position('right')
    axsA[1, 1].set_xscale(r'log')
    axsA[1, 1].set_xlim((1e0, 1e3))
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figA.tight_layout()
    figA.savefig('spe11b_time_series_boxA.png', bbox_inches='tight', dpi=300)

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
    axsB[0, 0].set_ylabel(r'mass [kt]')
    axsB[0, 0].set_xticklabels([])
    axsB[0, 0].set_xscale(r'log')
    axsB[0, 0].set_xlim((1e1, 1e3))
    axsB[0, 0].set_ylim((0, 4))
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
    axsB[0, 1].set_xticklabels([])
    axsB[0, 1].set_ylabel(r'mass [kt]')
    axsB[0, 1].yaxis.tick_right()
    axsB[0, 1].yaxis.set_label_position('right')
    axsB[0, 1].set_xscale(r'log')
    axsB[0, 1].set_xlim((1e1, 1e3))
    axsB[0, 1].set_ylim((0, 3))
    axsB[1, 0].set_title(r'Box B: dissolved CO2')
    axsB[1, 0].set_xlabel(r'time [y]')
    axsB[1, 0].set_ylabel(r'mass [kt]')
    axsB[1, 0].set_xscale(r'log')
    axsB[1, 0].set_xlim((1e1, 1e3))
    axsB[1, 0].set_ylim((0, 3))
    axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
    axsB[1, 1].set_xlabel(r'time [y]')
    axsB[1, 1].set_ylabel(r'mass [kt]')
    axsB[1, 1].yaxis.tick_right()
    axsB[1, 1].yaxis.set_label_position('right')
    axsB[1, 1].set_xscale(r'log')
    axsB[1, 1].set_xlim((1e1, 1e3))
    axsB[1, 1].set_ylim((0, 0.1))
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figB.tight_layout()
    figB.savefig('spe11b_time_series_boxB.png', bbox_inches='tight', dpi=300)

    axsC.set_title(r'Box C: convection')
    axsC.set_xlabel(r'time [y]')
    axsC.set_ylabel(r'$M$ [km]')
    axsC.set_xscale(r'log')
    axsC.set_xlim((1e1, 1e3))
    axsC.set_ylim((0, 20))
    axsC.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    figC.savefig('spe11b_time_series_boxC.png', bbox_inches='tight', dpi=300)

    axsT[0].set_title(r'CO2 in sealing units')
    axsT[0].set_xlabel(r'time [y]')
    axsT[0].set_ylabel(r'mass [t]')
    axsT[0].set_ylim(0, 5e2)
    axsT[0].set_xscale(r'log')
    axsT[0].set_xlim((1e0, 1e3))
    axsT[1].set_title(r'CO2 in boundary volumes')
    axsT[1].set_xlabel(r'time [y]')
    axsT[1].set_ylabel(r'mass [t]')
    axsT[1].set_ylim(0, 1e2)
    axsT[1].set_xscale(r'log')
    axsT[1].yaxis.tick_right()
    axsT[1].yaxis.set_label_position('right')
    axsT[1].set_xlim((4e2, 1e3))
    handles, labels = axsT[1].get_legend_handles_labels()
    figT.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figT.tight_layout()
    figT.savefig('spe11b_time_series_seal.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    assembleTimeSeries()
