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

def calculateConvection():
    """Calculate the convection integral for Case A of the 11th SPE CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This calculates the convection integral ((17) in the description) "
                    "based on given spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    fig, axs = plt.subplots(figsize=(5, 3))

    xSpace = np.arange(0.0, 2.8 + 5.0e-3, 1.0e-2)
    ySpace = np.arange(0.0, 1.2 + 5.0e-3, 1.0e-2)
    x, y = np.meshgrid(xSpace, ySpace)
    nX = xSpace.size - 1
    nY = ySpace.size - 1
    deltaX = deltaY = 1.0e-2

    for i, group in zip(range(len(groups)), groups):
        color = f'C{i}'

        if groupFolders:
            baseFolder = groupFolders[i]

        integral = []
        for h in range(121):
            if group[-2] != '-':
                if not groupFolders:
                    baseFolder = os.path.join(folder, group.lower())
                if group.lower() in groups_and_colors:
                    color = groups_and_colors[group.lower()]
                ls = '-'
            else:
                if not groupFolders:
                    baseFolder = os.path.join(folder, group[:-2].lower(), f'result{group[-1]}')
                if group[:-2].lower() in groups_and_colors:
                    color = groups_and_colors[group[:-2].lower()]
                if group[-1] == '1': ls = '-'
                elif group[-1] == '2': ls = '--'
                elif group[-1] == '3': ls = '-.'
                elif group[-1] == '4': ls = ':'

            fileName = os.path.join(baseFolder, f'spe11a_spatial_map_{h}h.csv')
            p, s, mCO2, mH2O, rhoG, rhoL, tmCO2 = getFieldValues(fileName, nX, nY)
            mCO2InBoxC = mCO2[9:41, 109:261]
            pInBoxC = p[9:41, 109:261]
            nXBoxC = len(mCO2InBoxC[0])
            nYBoxC = len(mCO2InBoxC)

            # Calculate CO2 solubility for the mean pressure in Box C at T = 20 °C.
            pInBoxC = np.nan_to_num(pInBoxC)
            pMean = np.mean(pInBoxC)
            A = solubility.computeA(293.15, pMean)
            B = solubility.computeB(293.15, pMean)
            y_H2O = (1 - B)/(1/A - B)
            xCO2_mol_mol = B*(1 - y_H2O)
            xCO2_kg_kg = 44/18*xCO2_mol_mol # convert from mol/mol to kg/kg

            gradX = 0.5/deltaX/xCO2_kg_kg*(mCO2InBoxC[1:nYBoxC-1, 2:nXBoxC] - mCO2InBoxC[1:nYBoxC-1, 0:nXBoxC-2])
            gradY = 0.5/deltaY/xCO2_kg_kg*(mCO2InBoxC[2:nYBoxC, 1:nXBoxC-1] - mCO2InBoxC[0:nYBoxC-2, 1:nXBoxC-1])
            gradX = np.nan_to_num(gradX)
            gradY = np.nan_to_num(gradY)
            norm = np.sqrt((np.square(gradX) + np.square(gradY)))
            integral.append(deltaX*deltaY*np.sum(norm))

        axs.plot(range(121), integral, label=group, color=color, linestyle=ls)

    axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs.set_title(r'Box C: convection from spatial maps')
    axs.set_xlabel(r'time [h]')
    axs.set_ylabel(r'$M$ [m]')
    axs.set_xlim(1e0, 7260.0/60)
    axs.set_ylim((0, 3.5))
    axs.set_xscale('log')
    fig.savefig('spe11a_convection_from_spatial_maps.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    calculateConvection()
