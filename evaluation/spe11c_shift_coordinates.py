# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3

""""
Script to shift x,y coordinates by half a cell in positive coordinate directions
"""

import os
import sys
import argparse
import numpy as np
np.set_printoptions(threshold=sys.maxsize)

def shiftCoordinates(fileName, newFolder):
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

    dx = dy = 25
    csvData[:, 0] = csvData[:, 0] + dx
    csvData[:, 1] = csvData[:, 1] + dy

    newFileName = os.path.join(newFolder, os.path.basename(fileName))
    header = 'x [m],y [m],z [m],WATER_PRESSURE [Pa],gas saturation [-],mass fraction of CO2 in liquid [-],mass fraction of H20 in vapor [-],phase mass density gas [kg/m3] ,phase mass density water [kg/m3],total mass CO2 [kg],temperature [Celsius]'
    np.savetxt(newFileName, csvData, delimiter=",", fmt='%.6e', header=header)

def turnZCoordinates():
    """Shift x,y coordinates by half a cell in positive coordinate directions for Case C of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script shifts x,y coordinates by half a cell in positive coordinate directions."
    )

    parser.add_argument('-f','--folder', help='path to folder containing spatial maps')

    cmdArgs = vars(parser.parse_args())
    folder = cmdArgs["folder"]

    newFolder = os.path.join(folder, 'shifted') 
    if not os.path.exists(newFolder):
        os.makedirs(newFolder)

    times = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
    for time in times:
        fileName = os.path.join(folder, f'spe11c_spatial_map_{time}y.csv')
        shiftCoordinates(fileName, newFolder)


if __name__ == "__main__":
    turnZCoordinates()
