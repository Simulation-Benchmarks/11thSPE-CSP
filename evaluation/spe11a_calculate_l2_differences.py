# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib

def myStr(value):
    return f'{value:.4e}'

def getFieldValues(fileName, nX, nY):
    p = np.zeros(nX*nY); p[:] = np.nan
    tmCO2 = np.zeros(nX*nY); tmCO2[:] = np.nan

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning nans.')
        return p, tmCO2

    skip_header = 0
    with open(fileName, "r") as file:
        if not (file.readline()[0]).isnumeric():
            skip_header = 1

    delimiter = ','

    csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=5)
    ind = np.lexsort((csvData[:,0], csvData[:,1]))
    csvData = csvData[ind]

    p[:] = csvData[:, 2]
    p[p < 1e0] = np.nan
    tmCO2[:] = csvData[:, 8]

    return p, tmCO2

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
    cellVolume = deltaX*deltaY

    # select file that contains impermeable cells with 'nan' pressure values
    fileNameSLB = os.path.join(folder, 'slb', 'spe11a', 'result1', f'spe11a_spatial_map_{hour}h.csv')
    pSLB, tmCO2 = getFieldValues(fileNameSLB, nX, nY)

    p = np.zeros([nX*nY, numGroups])
    tmCO2 = np.zeros([nX*nY, numGroups])

    for i, group in zip(range(numGroups), groups):
        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group, 'spe11a')
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1], 'spe11a', f'result{group[-1]}')

        fileName = os.path.join(baseFolder, f'spe11a_spatial_map_{hour}h.csv')
        p[:, i], tmCO2[:, i] = getFieldValues(fileName, nX, nY)

        # set values to 'nan' for impermeable cells
        p[:, i][np.isnan(pSLB)] = np.nan
        tmCO2[:, i][np.isnan(pSLB)] = np.nan

    l2NormP = np.zeros((numGroups, numGroups))
    l2SemiNormP = np.zeros((numGroups, numGroups))
    l2NormM = np.zeros((numGroups, numGroups))

    for i in range(numGroups-1):
        pI = p[:, i]
        tmCO2I = tmCO2[:, i]

        pIMean = np.nanmean(pI)
        pIVariation = pI - pIMean

        for j in range(i+1, numGroups):
            pJ = p[:, j]
            tmCO2J = tmCO2[:, j]

            pDiff = pI - pJ
            # set difference to zero for impermeable cells
            pDiff = np.nan_to_num(pDiff)
            l2NormP[i, j] = l2NormP[j, i] = np.sqrt(cellVolume*np.sum(np.square(pDiff)))

            tmCO2Diff = tmCO2I - tmCO2J
            tmCO2Diff = np.nan_to_num(tmCO2Diff)
            l2NormM[i, j] = l2NormM[j, i] = np.sqrt(cellVolume*np.sum(np.square(tmCO2Diff)))

            pJMean = np.nanmean(pJ)
            pJVariation = pJ - pJMean
            pVariationDiff = pIVariation - pJVariation
            pVariationDiff = np.nan_to_num(pVariationDiff)
            l2SemiNormP[i, j] = l2SemiNormP[j, i] = np.sqrt(cellVolume*np.sum(np.square(pVariationDiff)))

    with open(f'spe11a_pressure_l2_diff_{hour}h.csv', 'w') as f:
        print('#, ', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2NormP[i])), file=f)

    with open(f'spe11a_pressure_l2semi_diff_{hour}h.csv', 'w') as f:
        print('#, ', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2SemiNormP[i])), file=f)

    with open(f'spe11a_tmco2_l2_diff_{hour}h.csv', 'w') as f:
        print('#, ', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2NormM[i])), file=f)


if __name__ == "__main__":
    calculateL2Differences()
