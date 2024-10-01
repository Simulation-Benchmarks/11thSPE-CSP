# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
from spe11b_visualize_spatial_maps import getFieldValues

def calculateL2Differences():
    """Calculate the L2 pressure differences for Case B of the 11th SPE CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This calculates the L2 pressure differences based on given spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-t','--time', help='time in years', required=True)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    year = cmdArgs["time"]

    numGroups = len(groups)
    nX = 840
    nY = 120
    deltaX = deltaY = 1.0e1

    # select file that contains impermeable cells with 'nan' pressure values
    fileNameSLB = os.path.join(folder, 'slb', 'spe11b', f'spe11b_spatial_map_{year}y.csv')
    pSLB, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileNameSLB, nX, nY)

    l2Norm = np.zeros((numGroups, numGroups))
    h1SemiNorm = np.zeros((numGroups, numGroups))

    for i, groupI in zip(range(numGroups), groups):
        if groupFolders:
            baseFolderI = groupFolders[i]

        if groupI[-2] != '-':
            if not groupFolders:
                baseFolderI = os.path.join(folder, groupI.lower(), 'spe11b')
        else:
            if not groupFolders:
                baseFolderI = os.path.join(folder, groupI[:-2].lower(), 'spe11b', f'result{groupI[-1]}')

        fileNameI = os.path.join(baseFolderI, f'spe11b_spatial_map_{year}y.csv')
        pI, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileNameI, nX, nY)

        # set values to 'nan' for impermeable cells
        pI[pSLB == float('nan')] = float('nan')

        gradXI = 0.5/deltaX*(pI[1:nY-1, 2:nX] - pI[1:nY-1, 0:nX-2])
        gradYI = 0.5/deltaY*(pI[2:nY, 1:nX-1] - pI[0:nY-2, 1:nX-1])

        for j, groupJ in zip(range(numGroups), groups):
            if j <= i:
                continue

            if groupFolders:
                baseFolderJ = groupFolders[j]

            if groupJ[-2] != '-':
                if not groupFolders:
                    baseFolderJ = os.path.join(folder, groupJ.lower(), 'spe11b')
            else:
                if not groupFolders:
                    baseFolderJ = os.path.join(folder, groupJ[:-2].lower(), 'spe11b', f'result{groupJ[-1]}')

            fileNameJ = os.path.join(baseFolderJ, f'spe11b_spatial_map_{year}y.csv')
            pJ, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileNameJ, nX, nY)

            pDiff = pI - pJ
            # set difference to zero for impermeable cells
            pDiff = np.nan_to_num(pDiff)
            l2Norm[i, j] = l2Norm[j, i] = np.sqrt(deltaX*deltaY*np.sum(np.square(pDiff)))

            gradXJ = 0.5/deltaX*(pJ[1:nY-1, 2:nX] - pJ[1:nY-1, 0:nX-2])
            gradYJ = 0.5/deltaY*(pJ[2:nY, 1:nX-1] - pJ[0:nY-2, 1:nX-1])
            gradXDiff = gradXI - gradXJ
            gradYDiff = gradYI - gradYJ
            # set difference to zero for impermeable and their neighboring cells
            gradXDiff = np.nan_to_num(gradXDiff)
            gradYDiff = np.nan_to_num(gradYDiff)
            h1SemiNorm[i, j] = h1SemiNorm[j, i] = np.sqrt(deltaX*deltaY*(np.sum(np.square(gradXDiff)) + np.sum(np.square(gradYDiff))))


    for i, groupI in zip(range(numGroups), groups):
        print(groupI, l2Norm[i])

    for i, groupI in zip(range(numGroups), groups):
        print(groupI, h1SemiNorm[i])


if __name__ == "__main__":
    calculateL2Differences()
