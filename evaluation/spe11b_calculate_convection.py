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

def calculateConvection():
    """Calculate the convection integral for Case B of the 11th SPE CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This calculates the convection integral ((17) in the description) "
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

    fig, axs = plt.subplots(figsize=(5, 3))

    numGroups = len(groups)
    numTimeSteps = cmdArgs["numtimesteps"]
    dt = cmdArgs["timestepsize"]
    secondsPerYear = 3600*24*365
    table = np.zeros((numTimeSteps+1, numGroups+1))
    table[:, 0] = range(0, (numTimeSteps*dt+1)*secondsPerYear, dt*secondsPerYear)

    nX = 840
    nY = 120
    deltaX = deltaY = 1.0e1

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

        integral = []
        for year in range(0, numTimeSteps*dt+1, dt):
            fileName = os.path.join(baseFolder, f'spe11b_spatial_map_{year}y.csv')
            if not os.path.isfile(fileName):
                integral.append(np.nan)
                continue

            p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileName, nX, nY)
            mCO2InBoxC = mCO2[9:41, 329:781]
            pInBoxC = p[9:41, 329:781]
            tInBoxC = temp[9:41, 329:781]
            nXBoxC = len(mCO2InBoxC[0])
            nYBoxC = len(mCO2InBoxC)

            # Calculate CO2 solubility for the mean pressure and temperature in Box C.
            pMean = np.nanmean(pInBoxC)
            tMean = np.nanmean(tInBoxC)
            # The solubility functions expect temperature in Kelvin.
            A = solubility.computeA(tMean + 273.15, pMean)
            B = solubility.computeB(tMean + 273.15, pMean)
            y_H2O = (1 - B)/(1/A - B)
            xCO2_mol_mol = B*(1 - y_H2O)
            xCO2_kg_kg = 44*xCO2_mol_mol/(44*xCO2_mol_mol + 18*(1 - xCO2_mol_mol)) # convert from mol/mol to kg/kg

            gradX = 0.5/deltaX/xCO2_kg_kg*(mCO2InBoxC[1:nYBoxC-1, 2:nXBoxC] - mCO2InBoxC[1:nYBoxC-1, 0:nXBoxC-2])
            gradY = 0.5/deltaY/xCO2_kg_kg*(mCO2InBoxC[2:nYBoxC, 1:nXBoxC-1] - mCO2InBoxC[0:nYBoxC-2, 1:nXBoxC-1])
            gradX = np.nan_to_num(gradX)
            gradY = np.nan_to_num(gradY)
            norm = np.sqrt((np.square(gradX) + np.square(gradY)))
            # scale integral values to km for plotting
            integral.append(1e-3*deltaX*deltaY*np.sum(norm))

        axs.plot(range(0, numTimeSteps*dt+1, dt), integral, label=group, color=color, linestyle=ls)

        table[:, i+1] = integral
        # scale to m for file reporting
        table[:, i+1] = 1e3*table[:, i+1]

    axs.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    axs.set_title(r'Box C: convection from spatial maps')
    axs.set_xlabel(r'time [y]')
    axs.set_ylabel(r'$M$ [km]')
    axs.set_xscale('log')
    fig.savefig('spe11b_time_series_boxC_from_spatial_maps.png', bbox_inches='tight', dpi=300)

    np.savetxt('spe11b_mC_from_spatial_maps.csv', table, fmt='%.5e', delimiter=', ', header=header)

if __name__ == "__main__":
    calculateConvection()
