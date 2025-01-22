# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
from spe11a_visualize_spatial_maps import getFieldValues

def myStr(value):
    return f'{value:.4e}'

def calculateL2Differences():
    """Calculate the L2 pressure differences for Case A of the 11th SPE CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This calculates the L2 pressure differences based on given spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-t','--time', help='time in hours', required=True)

    cmdArgs = vars(parser.parse_args())
    groups = [x.lower() for x in cmdArgs["groups"]]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    hour = cmdArgs["time"]

    numGroups = len(groups)
    nX = 280
    nY = 120
    deltaX = deltaY = 1.0e-2

    # select file that contains impermeable cells with 'nan' pressure values
    fileNameSLB = os.path.join(folder, 'slb', 'spe11a', 'result1', f'spe11a_spatial_map_{hour}h.csv')
    pSLB, s, mCO2, mH2O, rhoG, rhoL, tmCO2 = getFieldValues(fileNameSLB, nX, nY)

    l2NormP = np.zeros((numGroups, numGroups))
    l2SemiNormP = np.zeros((numGroups, numGroups))
    h1SemiNormP = np.zeros((numGroups, numGroups))
    l2NormM = np.zeros((numGroups, numGroups))

    for i, groupI in zip(range(numGroups), groups):
        if groupFolders:
            baseFolderI = groupFolders[i]

        if not groupI[-1].isnumeric():
            if not groupFolders:
                baseFolderI = os.path.join(folder, groupI, 'spe11a')
        else:
            if not groupFolders:
                baseFolderI = os.path.join(folder, groupI[:-1], 'spe11a', f'result{groupI[-1]}')

        fileNameI0 = os.path.join(baseFolderI, f'spe11a_spatial_map_0h.csv')
        pI, s, mCO2, mH2O, rhoG, rhoL, tmCO2I0 = getFieldValues(fileNameI0, nX, nY)

        fileNameI = os.path.join(baseFolderI, f'spe11a_spatial_map_{hour}h.csv')
        pI, s, mCO2, mH2O, rhoG, rhoL, tmCO2I = getFieldValues(fileNameI, nX, nY)

        # subtract possibly added initial artificial CO2 mass
        tmCO2I = tmCO2I - tmCO2I0
        tmCO2I[tmCO2I < 0] = 0

        # set values to 'nan' for impermeable cells
        pI[np.isnan(pSLB)] = np.nan
        tmCO2I[np.isnan(pSLB)] = np.nan

        pIMean = np.nanmean(pI)
        pIVariation = pI - pIMean

        gradXI = 0.5/deltaX*(pI[1:nY-1, 2:nX] - pI[1:nY-1, 0:nX-2])
        gradYI = 0.5/deltaY*(pI[2:nY, 1:nX-1] - pI[0:nY-2, 1:nX-1])

        for j, groupJ in zip(range(numGroups), groups):
            if j <= i:
                continue

            if groupFolders:
                baseFolderJ = groupFolders[j]

            if not groupJ[-1].isnumeric():
                if not groupFolders:
                    baseFolderJ = os.path.join(folder, groupJ, 'spe11a')
            else:
                if not groupFolders:
                    baseFolderJ = os.path.join(folder, groupJ[:-1], 'spe11a', f'result{groupJ[-1]}')

            fileNameJ0 = os.path.join(baseFolderJ, f'spe11a_spatial_map_0h.csv')
            pJ, s, mCO2, mH2O, rhoG, rhoL, tmCO2J0 = getFieldValues(fileNameJ0, nX, nY)

            fileNameJ = os.path.join(baseFolderJ, f'spe11a_spatial_map_{hour}h.csv')
            pJ, s, mCO2, mH2O, rhoG, rhoL, tmCO2J = getFieldValues(fileNameJ, nX, nY)

            # subtract possibly added initial artificial CO2 mass
            tmCO2J = tmCO2J - tmCO2J0
            tmCO2J[tmCO2J < 0] = 0

            pJ[np.isnan(pSLB)] = np.nan
            tmCO2J[np.isnan(pSLB)] = np.nan

            pDiff = pI - pJ
            # set difference to zero for impermeable cells
            pDiff = np.nan_to_num(pDiff)
            l2NormP[i, j] = l2NormP[j, i] = np.sqrt(deltaX*deltaY*np.sum(np.square(pDiff)))

            tmCO2Diff = tmCO2I - tmCO2J
            tmCO2Diff = np.nan_to_num(tmCO2Diff)
            l2NormM[i, j] = l2NormM[j, i] = np.sqrt(deltaX*deltaY*np.sum(np.square(tmCO2Diff)))

            pJMean = np.nanmean(pJ)
            pJVariation = pJ - pJMean
            pVariationDiff = pIVariation - pJVariation
            pVariationDiff = np.nan_to_num(pVariationDiff)
            l2SemiNormP[i, j] = l2SemiNormP[j, i] = np.sqrt(deltaX*deltaY*np.sum(np.square(pVariationDiff)))

            gradXJ = 0.5/deltaX*(pJ[1:nY-1, 2:nX] - pJ[1:nY-1, 0:nX-2])
            gradYJ = 0.5/deltaY*(pJ[2:nY, 1:nX-1] - pJ[0:nY-2, 1:nX-1])
            gradXDiff = gradXI - gradXJ
            gradYDiff = gradYI - gradYJ
            # set difference to zero for impermeable and their neighboring cells
            gradXDiff = np.nan_to_num(gradXDiff)
            gradYDiff = np.nan_to_num(gradYDiff)
            h1SemiNormP[i, j] = h1SemiNormP[j, i] = np.sqrt(deltaX*deltaY*(np.sum(np.square(gradXDiff)) + np.sum(np.square(gradYDiff))))

    with open(f'spe11a_pressure_l2_diff_{hour}h.csv', 'w') as f:
        print('#, ', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2NormP[i])), file=f)

    with open(f'spe11a_pressure_l2semi_diff_{hour}h.csv', 'w') as f:
        print('#, ', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2SemiNormP[i])), file=f)

    with open(f'spe11a_pressure_h1semi_diff_{hour}h.csv', 'w') as f:
        print('#, ', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, h1SemiNormP[i])), file=f)

    with open(f'spe11a_tmco2_l2_diff_{hour}h.csv', 'w') as f:
        print('#, ', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2NormM[i])), file=f)


if __name__ == "__main__":
    calculateL2Differences()
