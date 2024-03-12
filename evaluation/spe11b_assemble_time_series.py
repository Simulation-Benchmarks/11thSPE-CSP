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
    figPT, axsPT = plt.subplots(1, 2, figsize=(9, 3))
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
        if group == "Calgary":
            delimiter = ' '

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365
        if group == "IFPEN" or group == "OpenGoSim":
            t = t - 1000

        axsP[0].plot(t, csvData[:, 1]/1e5, label=group, color=color, linestyle=ls)
        axsP[1].plot(t, csvData[:, 2]/1e5, label=group, color=color, linestyle=ls)

        axsPT[0].plot(t, csvData[:, 1]/1e5, label=group, color=color, linestyle=ls)
        axsPT[1].plot(t, csvData[:, 2]/1e5, label=group, color=color, linestyle=ls)

        axsA[0, 0].plot(t, csvData[:, 3], label=group, color=color, linestyle=ls)
        axsA[0, 1].plot(t, csvData[:, 4], label=group, color=color, linestyle=ls)
        axsA[1, 0].plot(t, csvData[:, 5], label=group, color=color, linestyle=ls)
        if group != "Kiel":
            axsA[1, 1].plot(t, csvData[:, 6], label=group, color=color, linestyle=ls)

        if group == "Kiel":
            axsB[0, 0].plot(t, csvData[:, 6], label=group, color=color, linestyle=ls)
            axsB[0, 1].plot(t, csvData[:, 7], label=group, color=color, linestyle=ls)
            axsB[1, 0].plot(t, csvData[:, 8], label=group, color=color, linestyle=ls)
        else:
            axsB[0, 0].plot(t, csvData[:, 7], label=group, color=color, linestyle=ls)
            axsB[0, 1].plot(t, csvData[:, 8], label=group, color=color, linestyle=ls)
            axsB[1, 0].plot(t, csvData[:, 9], label=group, color=color, linestyle=ls)
            axsB[1, 1].plot(t, csvData[:, 10], label=group, color=color, linestyle=ls)

        if group == "Kiel":
            axsC.plot(t, csvData[:, 9], label=group, color=color, linestyle=ls)
        else:
            axsC.plot(t, csvData[:, 11], label=group, color=color, linestyle=ls)

        if group == "Kiel":
            axsT[0].plot(t, csvData[:, 10], label=group, color=color, linestyle=ls)
        else:
            axsT[0].plot(t, csvData[:, 12], label=group, color=color, linestyle=ls)

        if group == "Kiel":
            axsT[1].plot(t, csvData[:, 11], label=group, color=color, linestyle=ls)
        elif group != "GEOS" and group != "UT-GCCC":
            axsT[1].plot(t, csvData[:, 13], label=group, color=color, linestyle=ls)

    axsP[0].set_title(r'sensor 1')
    axsP[0].set_xlabel(r'time [y]')
    axsP[0].set_ylabel(r'pressure [bar]')
    axsP[1].set_title(r'sensor 2')
    axsP[1].set_xlabel(r'time [y]')
    handles, labels = axsP[1].get_legend_handles_labels()
    figP.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=5)
    figP.savefig('spe11b_time_series_pressure.png', bbox_inches='tight')

    axsPT[0].set_title(r'sensor 1')
    axsPT[0].set_xlabel(r'time [y]')
    axsPT[0].set_ylabel(r'pressure [bar]')
    axsPT[1].set_title(r'sensor 2')
    axsPT[1].set_xlabel(r'time [y]')
    handles, labels = axsPT[1].get_legend_handles_labels()
    figPT.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=5)
    figPT.savefig('spe11b_time_series_pressure_zoom_time.png', bbox_inches='tight')

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO2')
    axsA[0, 0].set_ylabel(r'mass [kg]')
    axsA[0, 0].set_ylim(0, 5e7)
    axsA[0, 0].set_xticklabels([])
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO2')
    axsA[0, 1].set_ylim(0, 6e6)
    axsA[0, 1].set_xticklabels([])
    axsA[1, 0].set_title(r'Box A: dissolved CO2')
    axsA[1, 0].set_xlabel(r'time [y]')
    axsA[1, 0].set_ylabel(r'mass [kg]')
    axsA[1, 0].set_ylim(0, 1e7)
    axsA[1, 1].set_title(r'Box A: CO2 in the seal facies')
    axsA[1, 1].set_xlabel(r'time [y]')
    axsA[1, 1].set_ylim(0, 2e5)
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5)
    figA.savefig('spe11b_time_series_boxA.png', bbox_inches='tight')

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
    axsB[0, 0].set_ylabel(r'mass [kg]')
    axsB[0, 0].set_xticklabels([])
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
    axsB[0, 1].set_xticklabels([])
    axsB[1, 0].set_title(r'Box B: dissolved CO2')
    axsB[1, 0].set_xlabel(r'time [y]')
    axsB[1, 0].set_ylabel(r'mass [kg]')
    axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
    axsB[1, 1].set_xlabel(r'time [y]')
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5)
    figB.savefig('spe11b_time_series_boxB.png', bbox_inches='tight')

    axsC.set_title(r'Box C: convection')
    axsC.set_xlabel(r'time [y]')
    axsC.set_ylabel(r'$M$ [m]')
    axsC.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    figC.savefig('spe11b_time_series_boxC.png', bbox_inches='tight')

    axsT[0].set_title(r'CO2 in sealing units')
    axsT[0].set_xlabel(r'time [y]')
    axsT[0].set_ylabel(r'mass [kg]')
    axsT[0].set_ylim(0, 5e5)
    axsT[1].set_title(r'CO2 in boundary volumes')
    axsT[1].set_xlabel(r'time [y]')
    axsT[1].set_ylabel(r'mass [kg]')
    axsT[1].set_ylim(0, 1e5)
    handles, labels = axsT[1].get_legend_handles_labels()
    figT.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=5)
    figT.savefig('spe11b_time_series_seal.png', bbox_inches='tight')

if __name__ == "__main__":
    assembleTimeSeries()
