# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3

""""
Script to fix reporting grid
"""

import os
import sys
import argparse
import numpy as np
np.set_printoptions(threshold=sys.maxsize)

def fixGrid(fileName, newFolder):
    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning 0 values.')

    skip_header = 0
    with open(fileName, "r") as file:
        if not (file.readline()[0]).isnumeric():
            skip_header = 1

    delimiter = ','

    csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
    csvData[:,1] = 1.2e3 - csvData[:,1]

    numX = 840
    numY = 120
    cellWidth = 10
    fixedData = np.zeros((numX*numY, 10))

    for j, y in zip(np.arange(0, numY), np.arange(0, cellWidth*numY)):
        top = csvData[j*numX:int((j+0.5)*numX), :]
        bottom = csvData[int((j+0.5)*numX):(j+1)*numX, :]
        average = 0.5*(top[:, 2:] + bottom[:, 2:])
        fixedData[j*numX:(j+1)*numX, 0] = np.arange(0.5*cellWidth, cellWidth*numX, cellWidth)
        fixedData[j*numX:(j+1)*numX, 1] = cellWidth*(numY - 0.5 - j)
        fixedData[j*numX:(j+1)*numX:2, 2:] = average
        fixedData[j*numX+1:(j+1)*numX:2, 2:] = average

    newFileName = os.path.join(newFolder, os.path.basename(fileName))
    header = 'x [m],z [m],WATER_PRESSURE [Pa],gas saturation [-],mass fraction of CO2 in liquid [-],mass fraction of H20 in vapor [-],phase mass density gas [kg/m3],phase mass density water [kg/m3],total mass CO2 [kg],temperature [Celsius]'
    np.savetxt(newFileName, fixedData, delimiter=",", fmt='%.6e', header=header)

def fixReportingGrid():
    """Fix reporting grid for Case B of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script fixes the reporting grid."
    )

    parser.add_argument('-f','--folder', help='path to folder containing spatial maps')

    cmdArgs = vars(parser.parse_args())
    folder = cmdArgs["folder"]

    newFolder = os.path.join(folder, 'fixed') 
    if not os.path.exists(newFolder):
        os.makedirs(newFolder)

    times = np.arange(0, 1001, 5)
    for time in times:
        fileName = os.path.join(folder, f'spe11b_spatial_map_{time}y.csv')
        fixGrid(fileName, newFolder)


if __name__ == "__main__":
    fixReportingGrid()
