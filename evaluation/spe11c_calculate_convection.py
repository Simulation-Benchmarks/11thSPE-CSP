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
import sys
sys.path.append('..')
import thermodynamics.make_solubility_table as solubility


def getFieldValues(fileName, nX, nY, nZ):
    p = np.zeros([nX, nY, nZ]); p[:] = np.nan
    s = np.zeros([nX, nY, nZ]); s[:] = np.nan
    mCO2 = np.zeros([nX, nY, nZ]); mCO2[:] = np.nan
    temp = np.zeros([nX, nY, nZ]); temp[:] = np.nan

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning nans.')
        return p, mCO2, temp

    skip_header = 0
    with open(fileName, "r") as file:
        if not (file.readline()[0]).isnumeric():
            skip_header = 1

    delimiter = ','

    try:
        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
    except:
        print(f'Reading into table from file {fileName} failed. Returning nans.')
        return p, mCO2, temp

    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=3)
    csvData[:,2] = np.around(csvData[:,2], decimals=5)
    ind = np.lexsort((csvData[:,2], csvData[:,1], csvData[:,0]))
    csvData = csvData[ind]

    if len(np.unique(csvData[:,0])) != nX or len(np.unique(csvData[:,1])) != nY:
        print(f'Cannot find unique x ({len(np.unique(csvData[:,0]))} instead of {nX})', end=' ')
        print(f'or y ({len(np.unique(csvData[:,1]))} instead of {nY}) coordinates. Returning nans.')
        return p, mCO2, temp

    p = np.reshape(csvData[:, 3], (nX, nY, nZ)) if len(csvData[0]) > 3 else np.nan
    s = np.reshape(csvData[:, 4], (nX, nY, nZ)) if len(csvData[0]) > 4 else np.nan
    mCO2 = np.reshape(csvData[:, 5], (nX, nY, nZ)) if len(csvData[0]) > 5 else np.nan
    temp = np.reshape(csvData[:, 10], (nX, nY, nZ)) if len(csvData[0]) > 10 else np.nan

    p[p < 1e0] = np.nan
    mCO2[s > 1 - 1e-5] = np.nan
    mCO2[np.isnan(s)] = np.nan
    return p, mCO2, temp


def calculateConvection():
    """Calculate the convection integral for Case C of the 11th SPE CSP"""

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
    groups = [x.lower() for x in cmdArgs["groups"]]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    fig, axs = plt.subplots(figsize=(5, 3))

    numGroups = len(groups)
    timeSteps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
    numTimeSteps = len(timeSteps)
    secondsPerYear = 3600*24*365
    table = np.zeros((numTimeSteps, numGroups+1))
    table[:, 0] = secondsPerYear*np.array(timeSteps)

    nX = 168
    nY = 100
    nZ = 120
    deltaX = deltaY = 50
    deltaZ = 10
    cellVolume = deltaX*deltaY*deltaZ

    header = 'time [s]'
    for i, group in zip(range(numGroups), groups):
        color = f'C{i}'
        header = header + ', ' + group

        if groupFolders:
            baseFolder = groupFolders[i]

        if group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group)
            if group[:-1] in groups_and_colors:
                color = groups_and_colors[group[:-1]]
            if group[-1] == '1': ls = '-'
            elif group[-1] == '2': ls = '--'
            elif group[-1] == '3': ls = '-.'
            elif group[-1] == '4': ls = ':'
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, f"{group}1")
            if group in groups_and_colors:
                color = groups_and_colors[group]
            ls = '-'

        integral = []
        for year in timeSteps:
            fileName = os.path.join(baseFolder, f'spe11c_spatial_map_{year}y.csv')
            if not os.path.isfile(fileName):
                integral.append(np.nan)
                continue

            p, mCO2, temp = getFieldValues(fileName, nX, nY, nZ)
            if np.isnan(p).all():
                integral.append(np.nan)
                continue
                
            mCO2InBoxC = mCO2[66:156, :, 13:43]
            pInBoxC = p[66:156, :, 13:43]
            tInBoxC = temp[66:156, :, 13:43]
            nXBoxC = len(mCO2InBoxC)
            nYBoxC = len(mCO2InBoxC[0])
            nZBoxC = len(mCO2InBoxC[0][0])

            # Calculate CO2 solubility for the mean pressure and temperature in Box C.
            pMean = np.nanmean(pInBoxC)
            tMean = np.nanmean(tInBoxC)
            # The solubility functions expect temperature in Kelvin.
            A = solubility.computeA(tMean + 273.15, pMean)
            B = solubility.computeB(tMean + 273.15, pMean)
            y_H2O = (1 - B)/(1/A - B)
            xCO2_mol_mol = B*(1 - y_H2O)
            xCO2_kg_kg = 44*xCO2_mol_mol/(44*xCO2_mol_mol + 18*(1 - xCO2_mol_mol)) # convert from mol/mol to kg/kg

            gradX = 0.5/deltaX/xCO2_kg_kg*(mCO2InBoxC[2:nXBoxC, 1:nYBoxC-1, 1:nZBoxC-1] - mCO2InBoxC[0:nXBoxC-2, 1:nYBoxC-1, 1:nZBoxC-1])
            gradY = 0.5/deltaY/xCO2_kg_kg*(mCO2InBoxC[1:nXBoxC-1, 2:nYBoxC, 1:nZBoxC-1] - mCO2InBoxC[1:nXBoxC-1, 0:nYBoxC-2, 1:nZBoxC-1])
            gradZ = 0.5/deltaZ/xCO2_kg_kg*(mCO2InBoxC[1:nXBoxC-1, 1:nYBoxC-1, 2:nZBoxC] - mCO2InBoxC[1:nXBoxC-1, 1:nYBoxC-1, 0:nZBoxC-2])
            gradX = np.nan_to_num(gradX)
            gradY = np.nan_to_num(gradY)
            gradZ = np.nan_to_num(gradZ)
            norm = np.sqrt((np.square(gradX) + np.square(gradY) + np.square(gradZ)))
            # scale integral values to km2 for plotting
            integral.append(1e-6*cellVolume*np.sum(norm))

        axs.plot(timeSteps, integral, label=group, color=color, linestyle=ls)
        axs.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
        axs.set_title(r'Box C: convection from spatial maps')
        axs.set_xlabel(r'time [y]')
        axs.set_ylabel(r'$M$ [km$^2$]')
        axs.set_xscale('log')
        fig.savefig('spe11c_time_series_boxC_from_spatial_maps.png', bbox_inches='tight', dpi=300)

        table[:, i+1] = integral
        # scale to m2 for file reporting
        table[:, i+1] = 1e6*table[:, i+1]
        np.savetxt('spe11c_mC_from_spatial_maps.csv', table, fmt='%.5e', delimiter=', ', header=header)


if __name__ == "__main__":
    calculateConvection()
