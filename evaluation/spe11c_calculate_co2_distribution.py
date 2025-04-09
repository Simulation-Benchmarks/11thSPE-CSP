# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
from spe11c_visualize_spatial_maps import getFieldValues
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors
import sys
sys.path.append('..')
import thermodynamics.make_solubility_table as solubility

def calculateCO2Distribution():
    """Calculate the CO2 distribution among phases in Boxes A and B for Case C of the 11th SPE CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This calculates the CO2 distribution among phases in Boxes A and B "
                    "based on given spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-nt','--numtimesteps', help='number of time steps excluding the initial one', required=False, type=int, default=25)

    cmdArgs = vars(parser.parse_args())
    groups = [x.lower() for x in cmdArgs["groups"]]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))

    numGroups = len(groups)
    timeSteps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
    numTimeSteps = cmdArgs["numtimesteps"] + 1
    timeSteps = timeSteps[:numTimeSteps]
    secondsPerYear = 3600*24*365
    mobileATable = np.zeros((numTimeSteps, numGroups+1))
    mobileATable[:, 0] = secondsPerYear*np.array(timeSteps)
    immobileATable = np.zeros((numTimeSteps, numGroups+1))
    immobileATable[:, 0] = secondsPerYear*np.array(timeSteps)
    dissolvedATable = np.zeros((numTimeSteps, numGroups+1))
    dissolvedATable[:, 0] = secondsPerYear*np.array(timeSteps)
    sealATable = np.zeros((numTimeSteps, numGroups+1))
    sealATable[:, 0] = secondsPerYear*np.array(timeSteps)
    mobileBTable = np.zeros((numTimeSteps, numGroups+1))
    mobileBTable[:, 0] = secondsPerYear*np.array(timeSteps)
    immobileBTable = np.zeros((numTimeSteps, numGroups+1))
    immobileBTable[:, 0] = secondsPerYear*np.array(timeSteps)
    dissolvedBTable = np.zeros((numTimeSteps, numGroups+1))
    dissolvedBTable[:, 0] = secondsPerYear*np.array(timeSteps)
    sealBTable = np.zeros((numTimeSteps, numGroups+1))
    sealBTable[:, 0] = secondsPerYear*np.array(timeSteps)

    nX = 168
    nY = 100
    nZ = 120
    deltaX = deltaY = 50
    deltaZ = 10
    cellVolume = deltaX*deltaY*deltaZ

    facies = np.flip(np.load("spe11c_facies.npy"), 0)
    facies = np.transpose(facies, (1, 2, 0))
    porosity = np.zeros((nX, nY, nZ))
    porosity[facies == 1] = 0.1
    porosity[facies == 2] = 0.2
    porosity[facies == 3] = 0.2
    porosity[facies == 4] = 0.2
    porosity[facies == 5] = 0.25
    porosity[facies == 6] = 0.35
    # transform the interval [0:750] in physical z-coordinate to [0:650] in reference w-coordinate
    phiInBoxA = porosity[66:166, :, 0:65]
    faciesInBoxA = facies[66:166, :, 0:65]
    # transform the interval [750:1350] in physical z-coordinate to [650:1200] in reference w-coordinate
    phiInBoxB = porosity[2:66, :, 65:120]
    faciesInBoxB = facies[2:66, :, 65:120]


    header = 'time [s]'
    for i, group in zip(range(numGroups), groups):
        color = f'C{i}'
        header = header + ', ' + group

        if groupFolders:
            baseFolder = groupFolders[i]

        mobileA = []
        immobileA = []
        dissolvedA = []
        sealA = []
        mobileB = []
        immobileB = []
        dissolvedB = []
        sealB = []
        for year in timeSteps:
            if not group[-1].isnumeric():
                if not groupFolders:
                    baseFolder = os.path.join(folder, f"{group}1")
                if group in groups_and_colors:
                    color = groups_and_colors[group]
                ls = '-'
            else:
                if not groupFolders:
                    baseFolder = os.path.join(folder, group)
                if group[:-1] in groups_and_colors:
                    color = groups_and_colors[group[:-1]]
                if group[-1] == '1': ls = '-'
                elif group[-1] == '2': ls = '--'
                elif group[-1] == '3': ls = '-.'
                elif group[-1] == '4': ls = ':'

            fileName = os.path.join(baseFolder, f'spe11c_spatial_map_{year}y.csv')
            if not os.path.isfile(fileName):
                mobileA.append(np.nan)
                immobileA.append(np.nan)
                dissolvedA.append(np.nan)
                sealA.append(np.nan)
                mobileB.append(np.nan)
                immobileB.append(np.nan)
                dissolvedB.append(np.nan)
                sealB.append(np.nan)
                continue

            p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileName, nX, nY, nZ)
            if np.isnan(p).all():
                mobileA.append(np.nan)
                immobileA.append(np.nan)
                dissolvedA.append(np.nan)
                sealA.append(np.nan)
                mobileB.append(np.nan)
                immobileB.append(np.nan)
                dissolvedB.append(np.nan)
                sealB.append(np.nan)
                continue
                
            mCO2InBoxA = np.nan_to_num(mCO2[66:166, :, 0:65])
            mH2OInBoxA = np.nan_to_num(mH2O[66:166, :, 0:65])
            rhoGInBoxA = np.nan_to_num(rhoG[66:166, :, 0:65])
            rhoLInBoxA = np.nan_to_num(rhoL[66:166, :, 0:65])
            sInBoxA = np.nan_to_num(s[66:166, :, 0:65])
            nXBoxA = mCO2InBoxA.shape[0]
            nYBoxA = mCO2InBoxA.shape[1]
            nZBoxA = mCO2InBoxA.shape[2]
            mCO2InBoxB = np.nan_to_num(mCO2[2:66, :, 65:120])
            mH2OInBoxB = np.nan_to_num(mH2O[2:66, :, 65:120])
            rhoGInBoxB = np.nan_to_num(rhoG[2:66, :, 65:120])
            rhoLInBoxB = np.nan_to_num(rhoL[2:66, :, 65:120])
            sInBoxB = np.nan_to_num(s[2:66, :, 65:120])
            nXBoxB = mCO2InBoxB.shape[0]
            nYBoxB = mCO2InBoxB.shape[1]
            nZBoxB = mCO2InBoxB.shape[2]
            
            mCO2InGS = np.multiply(1 - mH2OInBoxA, sInBoxA)
            mCO2InGSRhoG = np.multiply(mCO2InGS, rhoGInBoxA)
            massGasInBoxA = np.multiply(mCO2InGSRhoG, phiInBoxA)

            mobMassInBoxA = np.zeros((nXBoxA, nYBoxA, nZBoxA))
            mobMassInBoxA[sInBoxA > 0.1] = massGasInBoxA[sInBoxA > 0.1]
            # integrate and scale to Mt for plotting
            mobileA.append(1e-9*cellVolume*np.sum(mobMassInBoxA))
            
            immMassInBoxA = np.zeros((nXBoxA, nYBoxA, nZBoxA))
            immMassInBoxA[sInBoxA <= 0.1] = massGasInBoxA[sInBoxA <= 0.1]
            # integrate and scale to Mt for plotting
            immobileA.append(1e-9*cellVolume*np.sum(immMassInBoxA))

            mCO2RhoL = np.multiply(mCO2InBoxA, rhoLInBoxA)
            mCO2RhoLS = np.multiply(mCO2RhoL, 1 - sInBoxA)
            dissMassInBoxA = np.multiply(mCO2RhoLS, phiInBoxA)
            # integrate and scale to Mt for plotting
            dissolvedA.append(1e-9*cellVolume*np.sum(dissMassInBoxA))

            sealMassInBoxA = np.zeros((nXBoxA, nYBoxA, nZBoxA))
            sealMassInBoxA[faciesInBoxA == 1] = massGasInBoxA[faciesInBoxA == 1] + dissMassInBoxA[faciesInBoxA == 1]
            # integrate and scale to Mt for plotting
            sealA.append(1e-9*cellVolume*np.sum(sealMassInBoxA))

            mCO2InGS = np.multiply(1 - mH2OInBoxB, sInBoxB)
            mCO2InGSRhoG = np.multiply(mCO2InGS, rhoGInBoxB)
            massGasInBoxB = np.multiply(mCO2InGSRhoG, phiInBoxB)

            mobMassInBoxB = np.zeros((nXBoxB, nYBoxB, nZBoxB))
            mobMassInBoxB[sInBoxB > 0.1] = massGasInBoxB[sInBoxB > 0.1]
            # integrate and scale to Mt for plotting
            mobileB.append(1e-9*cellVolume*np.sum(mobMassInBoxB))
            
            immMassInBoxB = np.zeros((nXBoxB, nYBoxB, nZBoxB))
            immMassInBoxB[sInBoxB <= 0.1] = massGasInBoxB[sInBoxB <= 0.1]
            # integrate and scale to Mt for plotting
            immobileB.append(1e-9*cellVolume*np.sum(immMassInBoxB))

            mCO2RhoL = np.multiply(mCO2InBoxB, rhoLInBoxB)
            mCO2RhoLS = np.multiply(mCO2RhoL, 1 - sInBoxB)
            dissMassInBoxB = np.multiply(mCO2RhoLS, phiInBoxB)
            # integrate and scale to Mt for plotting
            dissolvedB.append(1e-9*cellVolume*np.sum(dissMassInBoxB))

            sealMassInBoxB = np.zeros((nXBoxB, nYBoxB, nZBoxB))
            sealMassInBoxB[faciesInBoxB == 1] = massGasInBoxB[faciesInBoxB == 1] + dissMassInBoxB[faciesInBoxB == 1]
            # integrate and scale to Mt for plotting
            sealB.append(1e-9*cellVolume*np.sum(sealMassInBoxB))

        axsA[0,0].plot(timeSteps, mobileA, label=group, color=color, linestyle=ls)
        axsA[0,1].plot(timeSteps, immobileA, label=group, color=color, linestyle=ls)
        axsA[1,0].plot(timeSteps, dissolvedA, label=group, color=color, linestyle=ls)
        axsA[1,1].plot(timeSteps, sealA, label=group, color=color, linestyle=ls)
        axsB[0,0].plot(timeSteps, mobileB, label=group, color=color, linestyle=ls)
        axsB[0,1].plot(timeSteps, immobileB, label=group, color=color, linestyle=ls)
        axsB[1,0].plot(timeSteps, dissolvedB, label=group, color=color, linestyle=ls)
        axsB[1,1].plot(timeSteps, sealB, label=group, color=color, linestyle=ls)

        mobileATable[:, i+1] = mobileA
        immobileATable[:, i+1] = immobileA
        dissolvedATable[:, i+1] = dissolvedA
        sealATable[:, i+1] = sealA
        # scale to kg for file reporting
        mobileATable[:, i+1] = 1e9*mobileATable[:, i+1]
        immobileATable[:, i+1] = 1e9*immobileATable[:, i+1]
        dissolvedATable[:, i+1] = 1e9*dissolvedATable[:, i+1]
        sealATable[:, i+1] = 1e9*sealATable[:, i+1]
        mobileBTable[:, i+1] = mobileB
        immobileBTable[:, i+1] = immobileB
        dissolvedBTable[:, i+1] = dissolvedB
        sealBTable[:, i+1] = sealB
        # scale to kg for file reporting
        mobileBTable[:, i+1] = 1e9*mobileBTable[:, i+1]
        immobileBTable[:, i+1] = 1e9*immobileBTable[:, i+1]
        dissolvedBTable[:, i+1] = 1e9*dissolvedBTable[:, i+1]
        sealBTable[:, i+1] = 1e9*sealBTable[:, i+1]

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO$_2$')
    axsA[0, 0].set_ylabel(r'mass [Mt]')
    axsA[0, 0].set_xticklabels([])
    axsA[0, 0].set_xscale(r'log')
    axsA[0, 0].set_xlim((1e0, 1e3))
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO$_2$')
    axsA[0, 1].set_xticklabels([])
    axsA[0, 1].set_ylabel(r'mass [Mt]')
    axsA[0, 1].yaxis.tick_right()
    axsA[0, 1].yaxis.set_label_position('right')
    axsA[0, 1].set_xscale(r'log')
    axsA[0, 1].set_xlim((1e0, 1e3))
    axsA[1, 0].set_title(r'Box A: dissolved CO$_2$')
    axsA[1, 0].set_xlabel(r'time [y]')
    axsA[1, 0].set_ylabel(r'mass [Mt]')
    axsA[1, 0].set_xscale(r'log')
    axsA[1, 0].set_xlim((1e0, 1e3))
    axsA[1, 1].set_title(r'Box A: CO$_2$ in the seal facies')
    axsA[1, 1].set_xlabel(r'time [y]')
    axsA[1, 1].set_ylabel(r'mass [Mt]')
    axsA[1, 1].yaxis.tick_right()
    axsA[1, 1].yaxis.set_label_position('right')
    axsA[1, 1].set_xscale(r'log')
    axsA[1, 1].set_xlim((1e0, 1e3))
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figA.tight_layout()
    figA.savefig('spe11c_time_series_boxA_from_spatial_maps.png', bbox_inches='tight', dpi=300)

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO$_2$')
    axsB[0, 0].set_ylabel(r'mass [Mt]')
    axsB[0, 0].set_xticklabels([])
    axsB[0, 0].set_xscale(r'log')
    axsB[0, 0].set_xlim((1e0, 1e3))
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO$_2$')
    axsB[0, 1].set_xticklabels([])
    axsB[0, 1].set_ylabel(r'mass [Mt]')
    axsB[0, 1].yaxis.tick_right()
    axsB[0, 1].yaxis.set_label_position('right')
    axsB[0, 1].set_xscale(r'log')
    axsB[0, 1].set_xlim((1e0, 1e3))
    axsB[1, 0].set_title(r'Box B: dissolved CO$_2$')
    axsB[1, 0].set_xlabel(r'time [y]')
    axsB[1, 0].set_ylabel(r'mass [Mt]')
    axsB[1, 0].set_xscale(r'log')
    axsB[1, 0].set_xlim((1e0, 1e3))
    axsB[1, 1].set_title(r'Box B: CO$_2$ in the seal facies')
    axsB[1, 1].set_xlabel(r'time [y]')
    axsB[1, 1].set_ylabel(r'mass [Mt]')
    axsB[1, 1].yaxis.tick_right()
    axsB[1, 1].yaxis.set_label_position('right')
    axsB[1, 1].set_xscale(r'log')
    axsB[1, 1].set_xlim((1e0, 1e3))
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figB.tight_layout()
    figB.savefig('spe11c_time_series_boxB_from_spatial_maps.png', bbox_inches='tight', dpi=300)

    np.savetxt('spe11c_mobA_from_spatial_maps.csv', mobileATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11c_immA_from_spatial_maps.csv', immobileATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11c_dissA_from_spatial_maps.csv', dissolvedATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11c_sealA_from_spatial_maps.csv', sealATable, fmt='%.5e', delimiter=', ', header=header)

    np.savetxt('spe11c_mobB_from_spatial_maps.csv', mobileBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11c_immB_from_spatial_maps.csv', immobileBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11c_dissB_from_spatial_maps.csv', dissolvedBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11c_sealB_from_spatial_maps.csv', sealBTable, fmt='%.5e', delimiter=', ', header=header)

if __name__ == "__main__":
    calculateCO2Distribution()
