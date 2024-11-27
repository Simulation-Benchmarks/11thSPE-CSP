# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
from spe11a_visualize_spatial_maps import getFieldValues
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors
import sys
sys.path.append('..')
import thermodynamics.make_solubility_table as solubility

def calculateCO2Distribution():
    """Calculate the CO2 distribution among phases in Boxes A and B for Case A of the 11th SPE CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This calculates the CO2 distribution among phases in Boxes A and B "
                    "based on given spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-nt','--numtimesteps', help='number of time steps excluding the initial one', required=False, type=int, default=120)

    parser.add_argument('-dt','--timestepsize', help='size of the time step in hours between two spatial maps', required=False, type=int, default=1)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    figA, axsA = plt.subplots(2, 2, figsize=(9, 6))
    figB, axsB = plt.subplots(2, 2, figsize=(9, 6))

    numGroups = len(groups)
    numTimeSteps = cmdArgs["numtimesteps"]
    dt = cmdArgs["timestepsize"]
    secondsPerHour = 3600
    mobileATable = np.zeros((numTimeSteps+1, numGroups+1))
    mobileATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)
    immobileATable = np.zeros((numTimeSteps+1, numGroups+1))
    immobileATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)
    dissolvedATable = np.zeros((numTimeSteps+1, numGroups+1))
    dissolvedATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)
    sealATable = np.zeros((numTimeSteps+1, numGroups+1))
    sealATable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)
    mobileBTable = np.zeros((numTimeSteps+1, numGroups+1))
    mobileBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)
    immobileBTable = np.zeros((numTimeSteps+1, numGroups+1))
    immobileBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)
    dissolvedBTable = np.zeros((numTimeSteps+1, numGroups+1))
    dissolvedBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)
    sealBTable = np.zeros((numTimeSteps+1, numGroups+1))
    sealBTable[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerHour, dt*secondsPerHour)

    nX = 280
    nY = 120
    deltaX = deltaY = deltaZ = 1.0e-2

    facies = np.flipud(np.load("spe11b_facies.npy"))
    facies = facies[:, 1:nX*3:3]
    porosity = np.zeros((nY, nX))
    porosity[facies == 1] = 0.44
    porosity[facies == 2] = 0.43
    porosity[facies == 3] = 0.44
    porosity[facies == 4] = 0.45
    porosity[facies == 5] = 0.43
    porosity[facies == 6] = 0.46
    phiInBoxA = porosity[0:60, 109:280]
    faciesInBoxA = facies[0:60, 109:280]
    phiInBoxB = porosity[59:120, 0:110]
    faciesInBoxB = facies[59:120, 0:110]


    header = 'time [s]'
    for i, group in zip(range(numGroups), groups):
        if group[-2] == '-':
            header = header + ', ' + group[:-2].lower() + group[-1]
        else:
            header = header + ', ' + group.lower()
        color = f'C{i}'

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
        for hour in range(0, numTimeSteps*dt+1, dt):
            if group[-2] != '-':
                if not groupFolders:
                    baseFolder = os.path.join(folder, group.lower(), 'spe11a')
                if group.lower() in groups_and_colors:
                    color = groups_and_colors[group.lower()]
                ls = '-'
            else:
                if not groupFolders:
                    baseFolder = os.path.join(folder, group[:-2].lower(), 'spe11a', f'result{group[-1]}')
                if group[:-2].lower() in groups_and_colors:
                    color = groups_and_colors[group[:-2].lower()]
                if group[-1] == '1': ls = '-'
                elif group[-1] == '2': ls = '--'
                elif group[-1] == '3': ls = '-.'
                elif group[-1] == '4': ls = ':'

            fileName = os.path.join(baseFolder, f'spe11a_spatial_map_{hour}h.csv')
            if not os.path.isfile(fileName):
                mobileA.append(float('nan'))
                immobileA.append(float('nan'))
                dissolvedA.append(float('nan'))
                sealA.append(float('nan'))
                mobileB.append(float('nan'))
                immobileB.append(float('nan'))
                dissolvedB.append(float('nan'))
                sealB.append(float('nan'))
                continue

            p, s, mCO2, mH2O, rhoG, rhoL, tmCO2 = getFieldValues(fileName, nX, nY)
            mCO2InBoxA = np.nan_to_num(mCO2[0:60, 109:280])
            mH2OInBoxA = np.nan_to_num(mH2O[0:60, 109:280])
            rhoGInBoxA = np.nan_to_num(rhoG[0:60, 109:280])
            rhoLInBoxA = np.nan_to_num(rhoL[0:60, 109:280])
            sInBoxA = np.nan_to_num(s[0:60, 109:280])
            nXBoxA = len(mCO2InBoxA[0])
            nYBoxA = len(mCO2InBoxA)
            mCO2InBoxB = np.nan_to_num(mCO2[59:120, 0:110])
            mH2OInBoxB = np.nan_to_num(mH2O[59:120, 0:110])
            rhoGInBoxB = np.nan_to_num(rhoG[59:120, 0:110])
            rhoLInBoxB = np.nan_to_num(rhoL[59:120, 0:110])
            sInBoxB = np.nan_to_num(s[59:120, 0:110])
            nXBoxB = len(mCO2InBoxB[0])
            nYBoxB = len(mCO2InBoxB)
            
            mCO2InGS = np.multiply(1 - mH2OInBoxA, sInBoxA)
            mCO2InGSRhoG = np.multiply(mCO2InGS, rhoGInBoxA)
            massGasInBoxA = np.multiply(mCO2InGSRhoG, phiInBoxA)

            mobMassInBoxA = np.zeros((nYBoxA, nXBoxA))
            mobMassInBoxA[sInBoxA > 0.1] = massGasInBoxA[sInBoxA > 0.1]
            # integrate and scale to g for plotting
            mobileA.append(1e3*deltaX*deltaY*deltaZ*np.sum(mobMassInBoxA))
            
            immMassInBoxA = np.zeros((nYBoxA, nXBoxA))
            immMassInBoxA[sInBoxA <= 0.1] = massGasInBoxA[sInBoxA <= 0.1]
            # integrate and scale to g for plotting
            immobileA.append(1e3*deltaX*deltaY*deltaZ*np.sum(immMassInBoxA))

            mCO2RhoL = np.multiply(mCO2InBoxA, rhoLInBoxA)
            mCO2RhoLS = np.multiply(mCO2RhoL, 1 - sInBoxA)
            dissMassInBoxA = np.multiply(mCO2RhoLS, phiInBoxA)
            # integrate and scale to g for plotting
            dissolvedA.append(1e3*deltaX*deltaY*deltaZ*np.sum(dissMassInBoxA))

            sealMassInBoxA = np.zeros((nYBoxA, nXBoxA))
            sealMassInBoxA[faciesInBoxA == 1] = massGasInBoxA[faciesInBoxA == 1] + dissMassInBoxA[faciesInBoxA == 1]
            # integrate and scale to g for plotting
            sealA.append(1e3*deltaX*deltaY*deltaZ*np.sum(sealMassInBoxA))

            mCO2InGS = np.multiply(1 - mH2OInBoxB, sInBoxB)
            mCO2InGSRhoG = np.multiply(mCO2InGS, rhoGInBoxB)
            massGasInBoxB = np.multiply(mCO2InGSRhoG, phiInBoxB)

            mobMassInBoxB = np.zeros((nYBoxB, nXBoxB))
            mobMassInBoxB[sInBoxB > 0.1] = massGasInBoxB[sInBoxB > 0.1]
            # integrate and scale to g for plotting
            mobileB.append(1e3*deltaX*deltaY*deltaZ*np.sum(mobMassInBoxB))
            
            immMassInBoxB = np.zeros((nYBoxB, nXBoxB))
            immMassInBoxB[sInBoxB <= 0.1] = massGasInBoxB[sInBoxB <= 0.1]
            # integrate and scale to g for plotting
            immobileB.append(1e3*deltaX*deltaY*deltaZ*np.sum(immMassInBoxB))

            mCO2RhoL = np.multiply(mCO2InBoxB, rhoLInBoxB)
            mCO2RhoLS = np.multiply(mCO2RhoL, 1 - sInBoxB)
            dissMassInBoxB = np.multiply(mCO2RhoLS, phiInBoxB)
            # integrate and scale to g for plotting
            dissolvedB.append(1e3*deltaX*deltaY*deltaZ*np.sum(dissMassInBoxB))

            sealMassInBoxB = np.zeros((nYBoxB, nXBoxB))
            sealMassInBoxB[faciesInBoxB == 1] = massGasInBoxB[faciesInBoxB == 1] + dissMassInBoxB[faciesInBoxB == 1]
            # integrate and scale to g for plotting
            sealB.append(1e3*deltaX*deltaY*deltaZ*np.sum(sealMassInBoxB))

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
        mobileATable[:, i+1] = 1e-3*mobileATable[:, i+1]
        immobileATable[:, i+1] = 1e-3*immobileATable[:, i+1]
        dissolvedATable[:, i+1] = 1e-3*dissolvedATable[:, i+1]
        sealATable[:, i+1] = 1e-3*sealATable[:, i+1]
        mobileBTable[:, i+1] = mobileB
        immobileBTable[:, i+1] = immobileB
        dissolvedBTable[:, i+1] = dissolvedB
        sealBTable[:, i+1] = sealB
        # scale to kg for file reporting
        mobileBTable[:, i+1] = 1e-3*mobileBTable[:, i+1]
        immobileBTable[:, i+1] = 1e-3*immobileBTable[:, i+1]
        dissolvedBTable[:, i+1] = 1e-3*dissolvedBTable[:, i+1]
        sealBTable[:, i+1] = 1e-3*sealBTable[:, i+1]

    axsA[0, 0].set_title(r'Box A: mobile gaseous CO2')
    axsA[0, 0].set_ylabel(r'mass [g]')
    axsA[0, 0].set_xticklabels([])
    axsA[0, 0].set_xscale(r'log')
    axsA[0, 0].set_xlim((1e-1, 1.2e2))
    axsA[0, 1].set_title(r'Box A: immobile gaseous CO2')
    axsA[0, 1].set_xticklabels([])
    axsA[0, 1].set_ylabel(r'mass [g]')
    axsA[0, 1].yaxis.tick_right()
    axsA[0, 1].yaxis.set_label_position('right')
    axsA[0, 1].set_xscale(r'log')
    axsA[0, 1].set_xlim((1e-1, 1.2e2))
    axsA[1, 0].set_title(r'Box A: dissolved CO2')
    axsA[1, 0].set_xlabel(r'time [h]')
    axsA[1, 0].set_ylabel(r'mass [g]')
    axsA[1, 0].set_xscale(r'log')
    axsA[1, 0].set_xlim((1e-1, 1.2e2))
    axsA[1, 1].set_title(r'Box A: CO2 in the seal facies')
    axsA[1, 1].set_xlabel(r'time [h]')
    axsA[1, 1].set_ylabel(r'mass [g]')
    axsA[1, 1].yaxis.tick_right()
    axsA[1, 1].yaxis.set_label_position('right')
    axsA[1, 1].set_xscale(r'log')
    axsA[1, 1].set_xlim((1e-1, 1.2e2))
    handles, labels = axsA[1][1].get_legend_handles_labels()
    figA.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figA.tight_layout()
    figA.savefig('spe11a_time_series_boxA_from_spatial_maps.png', bbox_inches='tight', dpi=300)

    axsB[0, 0].set_title(r'Box B: mobile gaseous CO2')
    axsB[0, 0].set_ylabel(r'mass [g]')
    axsB[0, 0].set_xticklabels([])
    axsB[0, 0].set_xscale(r'log')
    axsB[0, 0].set_xlim((1e-1, 1.2e2))
    axsB[0, 1].set_title(r'Box B: immobile gaseous CO2')
    axsB[0, 1].set_xticklabels([])
    axsB[0, 1].set_ylabel(r'mass [g]')
    axsB[0, 1].yaxis.tick_right()
    axsB[0, 1].yaxis.set_label_position('right')
    axsB[0, 1].set_xscale(r'log')
    axsB[0, 1].set_xlim((1e-1, 1.2e2))
    axsB[1, 0].set_title(r'Box B: dissolved CO2')
    axsB[1, 0].set_xlabel(r'time [h]')
    axsB[1, 0].set_ylabel(r'mass [g]')
    axsB[1, 0].set_xscale(r'log')
    axsB[1, 0].set_xlim((1e-1, 1.2e2))
    axsB[1, 1].set_title(r'Box B: CO2 in the seal facies')
    axsB[1, 1].set_xlabel(r'time [h]')
    axsB[1, 1].set_ylabel(r'mass [g]')
    axsB[1, 1].yaxis.tick_right()
    axsB[1, 1].yaxis.set_label_position('right')
    axsB[1, 1].set_xscale(r'log')
    axsB[1, 1].set_xlim((1e-1, 1.2e2))
    handles, labels = axsB[1][1].get_legend_handles_labels()
    figB.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    figB.tight_layout()
    figB.savefig('spe11a_time_series_boxB_from_spatial_maps.png', bbox_inches='tight', dpi=300)

    np.savetxt('spe11a_mobile_boxA_from_spatial_maps.csv', mobileATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11a_immobile_boxA_from_spatial_maps.csv', immobileATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11a_dissolved_boxA_from_spatial_maps.csv', dissolvedATable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11a_seal_boxA_from_spatial_maps.csv', sealATable, fmt='%.5e', delimiter=', ', header=header)

    np.savetxt('spe11a_mobile_boxB_from_spatial_maps.csv', mobileBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11a_immobile_boxB_from_spatial_maps.csv', immobileBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11a_dissolved_boxB_from_spatial_maps.csv', dissolvedBTable, fmt='%.5e', delimiter=', ', header=header)
    np.savetxt('spe11a_seal_boxB_from_spatial_maps.csv', sealBTable, fmt='%.5e', delimiter=', ', header=header)

if __name__ == "__main__":
    calculateCO2Distribution()
