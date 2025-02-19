# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
from spe11b_visualize_spatial_maps import getFieldValues
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors
import sys
sys.path.append('..')
import thermodynamics.make_solubility_table as solubility

def calculateCO2Distribution():
    """Calculate the CO2 distribution among phases in Boxes A and B for Case B of the 11th SPE CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This calculates the CO2 distribution among phases in Boxes A and B "
                    "based on given spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-nt','--numtimesteps', help='number of time steps excluding the initial one', required=False, type=int, default=200)

    parser.add_argument('-dt','--timestepsize', help='size of the time step in years between two spatial maps', required=False, type=int, default=5)

    cmdArgs = vars(parser.parse_args())
    groups = [x.lower() for x in cmdArgs["groups"]]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))

    numGroups = len(groups)
    numTimeSteps = cmdArgs["numtimesteps"]
    dt = cmdArgs["timestepsize"]
    secondsPerYear = 3600*24*365
    mobileATable = np.zeros((numTimeSteps+1, numGroups+1))
    mobileATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)
    immobileATable = np.zeros((numTimeSteps+1, numGroups+1))
    immobileATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)
    dissolvedATable = np.zeros((numTimeSteps+1, numGroups+1))
    dissolvedATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)
    sealATable = np.zeros((numTimeSteps+1, numGroups+1))
    sealATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)
    mobileBTable = np.zeros((numTimeSteps+1, numGroups+1))
    mobileBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)
    immobileBTable = np.zeros((numTimeSteps+1, numGroups+1))
    immobileBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)
    dissolvedBTable = np.zeros((numTimeSteps+1, numGroups+1))
    dissolvedBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)
    sealBTable = np.zeros((numTimeSteps+1, numGroups+1))
    sealBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)

    nX = 840
    nY = 120
    deltaX = deltaY = 1.0e1

    facies = np.flipud(np.load("spe11b_facies.npy"))
    porosity = np.zeros((nY, nX))
    porosity[facies == 1] = 0.1
    porosity[facies == 2] = 0.2
    porosity[facies == 3] = 0.2
    porosity[facies == 4] = 0.2
    porosity[facies == 5] = 0.25
    porosity[facies == 6] = 0.35
    phiInBoxA = porosity[0:60, 329:830]
    faciesInBoxA = facies[0:60, 329:830]
    phiInBoxB = porosity[59:120, 99:330]
    faciesInBoxB = facies[59:120, 99:330]


    header = 'time [s]'
    for i, group in zip(range(numGroups), groups):
        color = f'C{i}'
        header = header + ', ' + group

        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group, 'spe11b')
            if group in groups_and_colors:
                color = groups_and_colors[group]
            ls = '-'
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1], 'spe11b', f'result{group[-1]}')
            if group[:-1] in groups_and_colors:
                color = groups_and_colors[group[:-1]]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'

        # get baseline of mass quantities
        fileName = os.path.join(baseFolder, f'spe11b_spatial_map_0y.csv')
        p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileName, nX, nY)
        mCO2InBoxA0 = np.nan_to_num(mCO2[0:60, 329:830])
        sInBoxA0 = np.nan_to_num(s[0:60, 329:830])
        nXBoxA = len(mCO2InBoxA0[0])
        nYBoxA = len(mCO2InBoxA0)
        mCO2InBoxB0 = np.nan_to_num(mCO2[59:120, 99:330])
        sInBoxB0 = np.nan_to_num(s[59:120, 99:330])
        nXBoxB = len(mCO2InBoxB0[0])
        nYBoxB = len(mCO2InBoxB0)

        mobileA = []
        immobileA = []
        dissolvedA = []
        sealA = []
        mobileB = []
        immobileB = []
        dissolvedB = []
        sealB = []
        for year in range(0, numTimeSteps*dt+1, dt):
            fileName = os.path.join(baseFolder, f'spe11b_spatial_map_{year}y.csv')
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

            p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileName, nX, nY)
            mCO2InBoxA = np.nan_to_num(mCO2[0:60, 329:830])
            mH2OInBoxA = np.nan_to_num(mH2O[0:60, 329:830]) - mCO2InBoxA0
            mCO2InBoxA[mCO2InBoxA < 0] = 0
            rhoGInBoxA = np.nan_to_num(rhoG[0:60, 329:830])
            rhoLInBoxA = np.nan_to_num(rhoL[0:60, 329:830])
            sInBoxA = np.nan_to_num(s[0:60, 329:830]) - sInBoxA0
            sInBoxA[sInBoxA < 0] = 0
            mCO2InBoxB = np.nan_to_num(mCO2[59:120, 99:330]) - mCO2InBoxB0
            mCO2InBoxB[mCO2InBoxB < 0] = 0
            mH2OInBoxB = np.nan_to_num(mH2O[59:120, 99:330])
            rhoGInBoxB = np.nan_to_num(rhoG[59:120, 99:330])
            rhoLInBoxB = np.nan_to_num(rhoL[59:120, 99:330])
            sInBoxB = np.nan_to_num(s[59:120, 99:330]) - sInBoxB0
            sInBoxB[sInBoxB < 0] = 0
            
            mCO2InGS = np.multiply(1 - mH2OInBoxA, sInBoxA)
            mCO2InGSRhoG = np.multiply(mCO2InGS, rhoGInBoxA)
            massGasInBoxA = np.multiply(mCO2InGSRhoG, phiInBoxA)

            mobMassInBoxA = np.zeros((nYBoxA, nXBoxA))
            mobMassInBoxA[sInBoxA > 0.1] = massGasInBoxA[sInBoxA > 0.1]
            # integrate and scale to kt for plotting
            mobileA.append(1e-6*deltaX*deltaY*np.sum(mobMassInBoxA))
            
            immMassInBoxA = np.zeros((nYBoxA, nXBoxA))
            immMassInBoxA[sInBoxA <= 0.1] = massGasInBoxA[sInBoxA <= 0.1]
            # integrate and scale to kt for plotting
            immobileA.append(1e-6*deltaX*deltaY*np.sum(immMassInBoxA))

            mCO2RhoL = np.multiply(mCO2InBoxA, rhoLInBoxA)
            mCO2RhoLS = np.multiply(mCO2RhoL, 1 - sInBoxA)
            dissMassInBoxA = np.multiply(mCO2RhoLS, phiInBoxA)
            # integrate and scale to kt for plotting
            dissolvedA.append(1e-6*deltaX*deltaY*np.sum(dissMassInBoxA))

            sealMassInBoxA = np.zeros((nYBoxA, nXBoxA))
            sealMassInBoxA[faciesInBoxA == 1] = massGasInBoxA[faciesInBoxA == 1] + dissMassInBoxA[faciesInBoxA == 1]
            # integrate and scale to kt for plotting
            sealA.append(1e-6*deltaX*deltaY*np.sum(sealMassInBoxA))

            mCO2InGS = np.multiply(1 - mH2OInBoxB, sInBoxB)
            mCO2InGSRhoG = np.multiply(mCO2InGS, rhoGInBoxB)
            massGasInBoxB = np.multiply(mCO2InGSRhoG, phiInBoxB)

            mobMassInBoxB = np.zeros((nYBoxB, nXBoxB))
            mobMassInBoxB[sInBoxB > 0.1] = massGasInBoxB[sInBoxB > 0.1]
            # integrate and scale to kt for plotting
            mobileB.append(1e-6*deltaX*deltaY*np.sum(mobMassInBoxB))
            
            immMassInBoxB = np.zeros((nYBoxB, nXBoxB))
            immMassInBoxB[sInBoxB <= 0.1] = massGasInBoxB[sInBoxB <= 0.1]
            # integrate and scale to kt for plotting
            immobileB.append(1e-6*deltaX*deltaY*np.sum(immMassInBoxB))

            mCO2RhoL = np.multiply(mCO2InBoxB, rhoLInBoxB)
            mCO2RhoLS = np.multiply(mCO2RhoL, 1 - sInBoxB)
            dissMassInBoxB = np.multiply(mCO2RhoLS, phiInBoxB)
            # integrate and scale to kt for plotting
            dissolvedB.append(1e-6*deltaX*deltaY*np.sum(dissMassInBoxB))

            sealMassInBoxB = np.zeros((nYBoxB, nXBoxB))
            sealMassInBoxB[faciesInBoxB == 1] = massGasInBoxB[faciesInBoxB == 1] + dissMassInBoxB[faciesInBoxB == 1]
            # integrate and scale to kt for plotting
            sealB.append(1e-6*deltaX*deltaY*np.sum(sealMassInBoxB))

        axsA[0,0].plot(range(0, numTimeSteps*dt+1, dt), mobileA, label=group, color=color, linestyle=ls)
        axsA[0,1].plot(range(0, numTimeSteps*dt+1, dt), immobileA, label=group, color=color, linestyle=ls)
        axsA[1,0].plot(range(0, numTimeSteps*dt+1, dt), dissolvedA, label=group, color=color, linestyle=ls)
        axsA[1,1].plot(range(0, numTimeSteps*dt+1, dt), sealA, label=group, color=color, linestyle=ls)
        axsB[0,0].plot(range(0, numTimeSteps*dt+1, dt), mobileB, label=group, color=color, linestyle=ls)
        axsB[0,1].plot(range(0, numTimeSteps*dt+1, dt), immobileB, label=group, color=color, linestyle=ls)
        axsB[1,0].plot(range(0, numTimeSteps*dt+1, dt), dissolvedB, label=group, color=color, linestyle=ls)
        axsB[1,1].plot(range(0, numTimeSteps*dt+1, dt), sealB, label=group, color=color, linestyle=ls)

        mobileATable[:, i+1] = mobileA
        immobileATable[:, i+1] = immobileA
        dissolvedATable[:, i+1] = dissolvedA
        sealATable[:, i+1] = sealA
        # scale to kg for file reporting
        mobileATable[:, i+1] = 1e6*mobileATable[:, i+1]
        immobileATable[:, i+1] = 1e6*immobileATable[:, i+1]
        dissolvedATable[:, i+1] = 1e6*dissolvedATable[:, i+1]
        sealATable[:, i+1] = 1e6*sealATable[:, i+1]
        mobileBTable[:, i+1] = mobileB
        immobileBTable[:, i+1] = immobileB
        dissolvedBTable[:, i+1] = dissolvedB
        sealBTable[:, i+1] = sealB
        # scale to kg for file reporting
        mobileBTable[:, i+1] = 1e6*mobileBTable[:, i+1]
        immobileBTable[:, i+1] = 1e6*immobileBTable[:, i+1]
        dissolvedBTable[:, i+1] = 1e6*dissolvedBTable[:, i+1]
        sealBTable[:, i+1] = 1e6*sealBTable[:, i+1]

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO$_2$')
    axsA[0, 0].set_ylabel(r'mass [kt]')
    axsA[0, 0].set_xticklabels([])
    axsA[0, 0].set_xscale(r'log')
    axsA[0, 0].set_xlim((1e0, 1e3))
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO$_2$')
    axsA[0, 1].set_xticklabels([])
    axsA[0, 1].set_ylabel(r'mass [kt]')
    axsA[0, 1].yaxis.tick_right()
    axsA[0, 1].yaxis.set_label_position('right')
    axsA[0, 1].set_xscale(r'log')
    axsA[0, 1].set_xlim((1e0, 1e3))
    axsA[1, 0].set_title(r'Box A: dissolved CO$_2$')
    axsA[1, 0].set_xlabel(r'time [y]')
    axsA[1, 0].set_ylabel(r'mass [kt]')
    axsA[1, 0].set_xscale(r'log')
    axsA[1, 0].set_xlim((1e0, 1e3))
    axsA[1, 1].set_title(r'Box A: CO$_2$ in the seal facies')
    axsA[1, 1].set_xlabel(r'time [y]')
    axsA[1, 1].set_ylabel(r'mass [kt]')
    axsA[1, 1].yaxis.tick_right()
    axsA[1, 1].yaxis.set_label_position('right')
    axsA[1, 1].set_xscale(r'log')
    axsA[1, 1].set_xlim((1e0, 1e3))
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figA.tight_layout()
    figA.savefig('spe11b_time_series_boxA_from_spatial_maps.png', bbox_inches='tight', dpi=300)

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO$_2$')
    axsB[0, 0].set_ylabel(r'mass [kt]')
    axsB[0, 0].set_xticklabels([])
    axsB[0, 0].set_xscale(r'log')
    axsB[0, 0].set_xlim((1e0, 1e3))
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO$_2$')
    axsB[0, 1].set_xticklabels([])
    axsB[0, 1].set_ylabel(r'mass [kt]')
    axsB[0, 1].yaxis.tick_right()
    axsB[0, 1].yaxis.set_label_position('right')
    axsB[0, 1].set_xscale(r'log')
    axsB[0, 1].set_xlim((1e0, 1e3))
    axsB[1, 0].set_title(r'Box B: dissolved CO$_2$')
    axsB[1, 0].set_xlabel(r'time [y]')
    axsB[1, 0].set_ylabel(r'mass [kt]')
    axsB[1, 0].set_xscale(r'log')
    axsB[1, 0].set_xlim((1e0, 1e3))
    axsB[1, 1].set_title(r'Box B: CO$_2$ in the seal facies')
    axsB[1, 1].set_xlabel(r'time [y]')
    axsB[1, 1].set_ylabel(r'mass [kt]')
    axsB[1, 1].yaxis.tick_right()
    axsB[1, 1].yaxis.set_label_position('right')
    axsB[1, 1].set_xscale(r'log')
    axsB[1, 1].set_xlim((1e0, 1e3))
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figB.tight_layout()
    figB.savefig('spe11b_time_series_boxB_from_spatial_maps.png', bbox_inches='tight', dpi=300)

    np.savetxt('spe11b_mobA_from_spatial_maps.csv', mobileATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11b_immA_from_spatial_maps.csv', immobileATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11b_dissA_from_spatial_maps.csv', dissolvedATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11b_sealA_from_spatial_maps.csv', sealATable, fmt='%.5e', delimiter=', ', header=header)

    np.savetxt('spe11b_mobB_from_spatial_maps.csv', mobileBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11b_immB_from_spatial_maps.csv', immobileBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11b_dissB_from_spatial_maps.csv', dissolvedBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11b_sealB_from_spatial_maps.csv', sealBTable, fmt='%.5e', delimiter=', ', header=header)

if __name__ == "__main__":
    calculateCO2Distribution()
