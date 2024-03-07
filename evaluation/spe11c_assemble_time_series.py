#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def assembleTimeSeries():
    """Visualize time series for the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script visualizes the time series quantities "
                    "as required by the CSP description."
    )

    fileNames = ["/media/bernd/bernd/spe11/geos/spe11c_time_series.csv",
                 "/media/bernd/bernd/spe11/ifpen/spe11c_time_series.csv",
                 "/media/bernd/bernd/spe11/opengosim/spe11c_time_series.csv"]
    groups = ["GEOS", "IFPEN", "OpenGoSim"]
    colors = ["C0", "C2", "C3"]

    font = {'size' : 14}
    matplotlib.rc('font', **font)
    #plt.rcParams.update({
    #    "text.usetex": True,
    #    "font.family": "monospace",
    #    "legend.columnspacing": 1.0,
    #    "legend.handlelength": 1.2
    #})

    figP, axsP = plt.subplots(1, 2, figsize=(9, 3))
    figPT, axsPT = plt.subplots(1, 2, figsize=(9, 3))
    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))
    figC, axsC = plt.subplots(figsize=(5, 3))
    figT, axsT = plt.subplots(1, 2, figsize=(9, 3))

    for fileName, group, color in zip(fileNames, groups, colors):
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
        t = csvData[:, 0]/60/60/24/365
        if group == "IFPEN" or group == "OpenGoSim":
            t = t - 1000

        axsP[0].plot(t, csvData[:, 1]/1e5, label=group, color=color)
        axsP[0].set_title(r'sensor 1')
        axsP[0].set_xlabel(r'time [y]')
        axsP[0].set_ylabel(r'pressure [bar]')
        #axsP[0].set_ylim(1.09e0, 1.15e0)
        #axsP[0].set_xlim(-1.0/60, 7260.0/60)
        #axsP[0].set_yticks([1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15])

        axsP[1].plot(t, csvData[:, 2]/1e5, label=group, color=color)
        axsP[1].set_title(r'sensor 2')
        axsP[1].set_xlabel(r'time [y]')
        #axsP[1].set_ylim(1.03e0, 1.09e0)
        #axsP[1].set_xlim(-1.0/60, 7260.0/60)
        #axsP[1].set_yticks([1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09])

        axsPT[0].plot(t, csvData[:, 1]/1e5, label=group, color=color)
        axsPT[0].set_title(r'sensor 1')
        axsPT[0].set_xlabel(r'time [y]')
        axsPT[0].set_ylabel(r'pressure [bar]')
        #axsPT[0].set_xlim(-0.1/60, 610.0/60)
        #axsPT[0].set_ylim(1.09e0, 1.15e0)
        #axsPT[0].set_yticks([1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15])

        axsPT[1].plot(t, csvData[:, 2]/1e5, label=group, color=color)
        axsPT[1].set_title(r'sensor 2')
        axsPT[1].set_xlabel(r'time [y]')
        #axsPT[1].set_xlim(-0.1/60, 610.0/60)
        #axsPT[1].set_ylim(1.03e0, 1.09e0)
        #axsPT[1].set_yticks([1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09])

        axsA[0, 0].plot(t, csvData[:, 3], label=group, color=color)
        axsA[0, 0].set_title(r'Box A: mobile gaseous CO2')
        axsA[0, 0].set_ylabel(r'mass [kg]')
        #axsA[0, 0].set_xlim(-1.0/60, 7260.0/60)
        #axsA[0, 0].set_ylim(0, 5e7)
        axsA[0, 0].set_xticklabels([])

        axsA[0, 1].plot(t, csvData[:, 4], label=group, color=color)
        axsA[0, 1].set_title(r'Box A: immobile gaseous CO2')
        #axsA[0, 1].set_xlim(-1.0/60, 7260.0/60)
        #axsA[0, 1].set_ylim(0, 6e6)
        axsA[0, 1].set_xticklabels([])

        axsA[1, 0].plot(t, csvData[:, 5], label=group, color=color)
        axsA[1, 0].set_title(r'Box A: dissolved CO2')
        axsA[1, 0].set_xlabel(r'time [y]')
        axsA[1, 0].set_ylabel(r'mass [kg]')
        #axsA[1, 0].set_xlim(-1.0/60, 7260.0/60)
        #axsA[1, 0].set_ylim(0, 1e7)

        if group != "Kiel":
            axsA[1, 1].plot(t, csvData[:, 6], label=group, color=color)
        axsA[1, 1].set_title(r'Box A: CO2 in the seal facies')
        axsA[1, 1].set_xlabel(r'time [y]')
        #axsA[1, 1].set_xlim(-1.0/60, 7260.0/60)
        #axsA[1, 1].set_ylim(0, 2e5)

        if group == "Kiel":
            axsB[0, 0].plot(t, csvData[:, 6], label=group, color=color)
        else:
            axsB[0, 0].plot(t, csvData[:, 7], label=group, color=color)
        axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
        axsB[0, 0].set_ylabel(r'mass [kg]')
        #axsB[0, 0].set_xlim(-1.0/60, 7260.0/60)
        #axsB[0, 0].set_ylim(-0.01, 0.6)
        axsB[0, 0].set_xticklabels([])

        if group == "Kiel":
            axsB[0, 1].plot(t, csvData[:, 7], label=group, color=color)
        else:
            axsB[0, 1].plot(t, csvData[:, 8], label=group, color=color)
        axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
        #axsB[0, 1].set_xlim(-1.0/60, 7260.0/60)
        #axsB[0, 1].set_ylim(-0.001, 0.07)
        axsB[0, 1].set_xticklabels([])

        if group == "Kiel":
            axsB[1, 0].plot(t, csvData[:, 8], label=group, color=color)
        else:
            axsB[1, 0].plot(t, csvData[:, 9], label=group, color=color)
        axsB[1, 0].set_title(r'Box B: dissolved CO2')
        axsB[1, 0].set_xlabel(r'time [y]')
        axsB[1, 0].set_ylabel(r'mass [kg]')
        #axsB[1, 0].set_xlim(-1.0/60, 7260.0/60)
        #axsB[1, 0].set_ylim(-0.01, 2.5)

        if group != "Kiel":
            axsB[1, 1].plot(t, csvData[:, 10], label=group, color=color)
        axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
        axsB[1, 1].set_xlabel(r'time [y]')
        #axsB[1, 1].set_xlim(-1.0/60, 7260.0/60)
        #axsB[1, 1].set_ylim(-0.01, 0.6)

        if group == "Kiel":
            axsC.plot(t, csvData[:, 9], label=group, color=color)
        else:
            axsC.plot(t, csvData[:, 11], label=group, color=color)
        axsC.set_title(r'Box C: convection')
        axsC.set_xlabel(r'time [y]')
        axsC.set_ylabel(r'$M$ [m]')
        #axsC.set_xlim(-1.0/60, 7260.0/60)
        axsC.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if group == "Kiel":
            axsT[0].plot(t, csvData[:, 10], label=group, color=color)
        else:
            axsT[0].plot(t, csvData[:, 12], label=group, color=color)
        axsT[0].set_title(r'CO2 in sealing units')
        axsT[0].set_xlabel(r'time [y]')
        axsT[0].set_ylabel(r'mass [kg]')
        #axsT[0].set_xlim(-1.0/60, 7260.0/60)
        #axsT[0].set_ylim(0, 5e5)

        if group == "Kiel":
            axsT[1].plot(t, csvData[:, 11], label=group, color=color)
        elif group != "GEOS" and group != "UT-GCCC":
            axsT[1].plot(t, csvData[:, 13], label=group, color=color)
        axsT[1].set_title(r'CO2 in boundary volumes')
        axsT[1].set_xlabel(r'time [y]')
        axsT[1].set_ylabel(r'mass [kg]')
        #axsT[1].set_xlim(-1.0/60, 7260.0/60)
        #axsT[1].set_ylim(0, 1e5)

    handles, labels = axsP[1].get_legend_handles_labels()
    figP.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=5)
    figP.savefig('spe11c_time_series_pressure.png', bbox_inches='tight')

    handles, labels = axsPT[1].get_legend_handles_labels()
    figPT.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=5)
    figPT.savefig('spe11c_time_series_pressure_zoom_time.png', bbox_inches='tight')

    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5)
    figA.savefig('spe11c_time_series_boxA.png', bbox_inches='tight')

    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=5)
    figB.savefig('spe11c_time_series_boxB.png', bbox_inches='tight')

    figC.savefig('spe11c_time_series_boxC.png', bbox_inches='tight')

    handles, labels = axsT[1].get_legend_handles_labels()
    figT.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.22), ncol=5)
    figT.savefig('spe11c_time_series_seal.png', bbox_inches='tight')

if __name__ == "__main__":
    assembleTimeSeries()
