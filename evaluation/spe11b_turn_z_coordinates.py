# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3

""""
Script to turn z coordinates around
"""

import os
import sys
import argparse
import numpy as np
np.set_printoptions(threshold=sys.maxsize)

def turnCoordinates(fileName, newFolder):
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

    newFileName = os.path.join(newFolder, os.path.basename(fileName))
    header = 'x [m],z [m],WATER_PRESSURE [Pa],gas saturation [-],mass fraction of CO2 in liquid [-],mass fraction of H20 in vapor [-],phase mass density gas [kg/m3] ,phase mass density water [kg/m3],total mass CO2 [kg],temperature [Celsius]'
    np.savetxt(newFileName, csvData, delimiter=",", fmt='%.6e', header=header)

def turnZCoordinates():
    """Turn z coordinates around for Case B of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script turns z coordinates around."
    )

    parser.add_argument('-f','--folder', help='path to folder containing spatial maps')

    cmdArgs = vars(parser.parse_args())
    folder = cmdArgs["folder"]

    newFolder = os.path.join(folder, 'sorted') 
    if not os.path.exists(newFolder):
        os.makedirs(newFolder)

    times = np.arange(0, 1001, 5).tolist()
    for time in times:
        fileName = os.path.join(folder, f'spe11b_spatial_map_{time}y.csv')
        turnCoordinates(fileName, newFolder)


if __name__ == "__main__":
    turnZCoordinates()
