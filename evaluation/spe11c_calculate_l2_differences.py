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

def getFieldValues(fileName, nX, nY, nZ):
    p = np.zeros(nX*nY*nZ); p[:] = np.nan
    tmCO2 = np.zeros(nX*nY*nZ); tmCO2[:] = np.nan
    temp = np.zeros(nX*nY*nZ); temp[:] = np.nan

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning nans.')
        return p, tmCO2, temp

    skip_header = 0
    with open(fileName, "r") as file:
        if not (file.readline()[0]).isnumeric():
            skip_header = 1

    delimiter = ','

    csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=3)
    csvData[:,2] = np.around(csvData[:,2], decimals=5)
    ind = np.lexsort((csvData[:,2], csvData[:,1], csvData[:,0]))
    csvData = csvData[ind]

    if len(np.unique(csvData[:,0])) != nX or len(np.unique(csvData[:,1])) != nY:
        print('Cannot find unique x or y coordinates. Returning nans.')
        return p, tmCO2, temp

    p[:] = csvData[:, 3]
    p[p < 1e0] = float('nan')
    tmCO2[:] = csvData[:, 9]
    temp[:] = csvData[:, 10]

    return p, tmCO2, temp

def calculateL2Differences():
    """Calculate the L2 pressure differences for Case C of the 11th SPE CSP"""

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
    groups = [x.lower() for x in cmdArgs["groups"]]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    year = cmdArgs["time"]

    numGroups = len(groups)
    nX = 168
    nY = 100
    nZ = 120
    deltaX = deltaY = 50
    deltaZ = 10
    cellVolume = deltaX*deltaY*deltaZ

    # select file that contains impermeable cells with 'nan' pressure values
    fileNameSLB = os.path.join(folder, 'slb', 'spe11c', f'spe11c_spatial_map_{year}y.csv')
    pSLB, tmCO2, temp = getFieldValues(fileNameSLB, nX, nY, nZ)


    p = np.zeros([nX*nY*nZ, numGroups])
    tmCO2 = np.zeros([nX*nY*nZ, numGroups])
    temp = np.zeros([nX*nY*nZ, numGroups])

    for i, group in zip(range(numGroups), groups):
        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group, 'spe11c')
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1], 'spe11c', f'result{group[-1]}')

        fileName = os.path.join(baseFolder, f'spe11c_spatial_map_{year}y.csv')
        p[:, i], tmCO2[:, i], temp[:, i] = getFieldValues(fileName, nX, nY, nZ)

        # set values to 'nan' for impermeable cells
        p[:, i][np.isnan(pSLB)] = float('nan')
        tmCO2[:, i][np.isnan(pSLB)] = float('nan')

    l2NormP = np.zeros((numGroups, numGroups))
    l2SemiNormP = np.zeros((numGroups, numGroups))
    l2NormT = np.zeros((numGroups, numGroups))
    l2SemiNormT = np.zeros((numGroups, numGroups))
    l2NormM = np.zeros((numGroups, numGroups))

    for i in range(numGroups-1):
        pI = p[:, i]
        tmCO2I = tmCO2[:, i]
        tempI = temp[:, i]

        pIMean = np.nanmean(pI)
        pIVariation = pI - pIMean
        tempIMean = np.nanmean(tempI)
        tempIVariation = tempI - tempIMean

        for j in range(i+1, numGroups):
            pJ = p[:, j]
            tmCO2J = tmCO2[:, j]
            tempJ = temp[:, j]

            pDiff = pI - pJ
            # set difference to zero for impermeable cells
            pDiff = np.nan_to_num(pDiff)
            l2NormP[i, j] = l2NormP[j, i] = np.sqrt(cellVolume*np.sum(np.square(pDiff)))

            tempDiff = tempI - tempJ
            tempDiff = np.nan_to_num(tempDiff)
            l2NormT[i, j] = l2NormT[j, i] = np.sqrt(cellVolume*np.sum(np.square(tempDiff)))

            tmCO2Diff = tmCO2I - tmCO2J
            tmCO2Diff = np.nan_to_num(tmCO2Diff)
            l2NormM[i, j] = l2NormM[j, i] = np.sqrt(cellVolume*np.sum(np.square(tmCO2Diff)))

            pJMean = np.nanmean(pJ)
            pJVariation = pJ - pJMean
            pVariationDiff = pIVariation - pJVariation
            pVariationDiff = np.nan_to_num(pVariationDiff)
            l2SemiNormP[i, j] = l2SemiNormP[j, i] = np.sqrt(cellVolume*np.sum(np.square(pVariationDiff)))

            tempJMean = np.nanmean(tempJ)
            tempJVariation = tempJ - tempJMean
            tempVariationDiff = tempIVariation - tempJVariation
            tempVariationDiff = np.nan_to_num(tempVariationDiff)
            l2SemiNormT[i, j] = l2SemiNormT[j, i] = np.sqrt(cellVolume*np.sum(np.square(tempVariationDiff)))

    with open(f'spe11c_pressure_l2_diff_{year}y.csv', 'w') as f:
        print('#', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2NormP[i])), file=f)

    with open(f'spe11c_pressure_l2semi_diff_{year}y.csv', 'w') as f:
        print('#', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2SemiNormP[i])), file=f)

    with open(f'spe11c_temperature_l2_diff_{year}y.csv', 'w') as f:
        print('#', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2NormT[i])), file=f)

    with open(f'spe11c_temperature_l2semi_diff_{year}y.csv', 'w') as f:
        print('#', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2SemiNormT[i])), file=f)

    with open(f'spe11c_tmco2_l2_diff_{year}y.csv', 'w') as f:
        print('#', ', '.join(map(str, groups)), file=f)
        for i, groupI in zip(range(numGroups), groups):
            print(groupI + ',', ', '.join(map(myStr, l2NormM[i])), file=f)


if __name__ == "__main__":
    calculateL2Differences()
