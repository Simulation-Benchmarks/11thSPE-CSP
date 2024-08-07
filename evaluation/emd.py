# SPDX-FileCopyrightText: 2024 Jakub W. Both, Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
# Test script for the Earth movers distance to measure differences in distributions in 2d.
import ot
import os
import sys
import argparse
import numpy as np
from PIL import Image
import time
np.set_printoptions(threshold=sys.maxsize)

def getFieldValues(fileName, nX, nY):
    x = np.zeros([nY, nX])
    y = np.zeros([nY, nX])
    tmCO2 = np.zeros([nY, nX])

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning 0 values.')
        return tmCO2

    skip_header = 0
    with open(fileName, "r") as file:
        if not (file.readline()[0]).isnumeric():
            skip_header = 1

    csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)
    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=5)
    ind = np.lexsort((csvData[:,0], csvData[:,1]))
    csvData = csvData[ind]
    idx = 8
    for i in np.arange(0, nY):
        x[i, :] = csvData[i*nX:(i+1)*nX, 0]
        y[i, :] = csvData[i*nX:(i+1)*nX, 1]
        tmCO2[i, :] = csvData[i*nX:(i+1)*nX, idx] if len(csvData[0]) > idx else 0

    tmCO2[np.isnan(tmCO2)] = 0

    return x, y, tmCO2

def calculateEMD(inFileName1, inFileName2, nX, nY, shrinkingFactor):
    # NOTE: ot has severe memory restrictions. Furthermore, computing the exact EMD has
    # O(n**3) complexity and can therefore be quite slow.

    # Get CO2 mass distributions on the reporting grid.
    x1, y1, co2Mass1 = getFieldValues(inFileName1, nX, nY)
    x2, y2, co2Mass2 = getFieldValues(inFileName2, nX, nY)

    # Shrink by the selected factor.
    x1 = x1[::shrinkingFactor, ::shrinkingFactor]
    y1 = y1[::shrinkingFactor, ::shrinkingFactor]
    co2Mass1 = co2Mass1[::shrinkingFactor, ::shrinkingFactor]
    x2 = x2[::shrinkingFactor, ::shrinkingFactor]
    y2 = y2[::shrinkingFactor, ::shrinkingFactor]
    co2Mass2 = co2Mass2[::shrinkingFactor, ::shrinkingFactor]

    # Flatten the arrays and normalize.
    # ot requires a and b to be compatible, i.e., that their sums are equal.
    x1 = x1.flatten()
    y1 = y1.flatten()
    co2Mass1 = co2Mass1.flatten()/np.sum(co2Mass1)
    x2 = x2.flatten()
    y2 = y2.flatten()
    co2Mass2 = co2Mass2.flatten()/np.sum(co2Mass2)

    # Cell centers of all cells - x and y coordinates.
    cc1 = np.vstack((x1, y1)).T
    cc2 = np.vstack((x2, y2)).T
        
    # Calculate the distance matrix and apply the EMD algorithm.
    M = ot.dist(cc1, cc2, metric="euclidean")
    dist_ot = ot.emd2(co2Mass1, co2Mass2, M, numItermax=500000)

    print(f'{inFileName1}, {inFileName2}: {dist_ot}')

    return dist_ot

def emd():
    """Calculate the Wasserstein distance between the CO2 mass distributions of two spatial maps"""

    parser = argparse.ArgumentParser(
        description="This script calculates the Wasserstein distance between the CO2 mass distributions of two spatial maps."
    )
    parser.add_argument("-in1", "--infilename1", help="The file name of the first spatial map.")
    parser.add_argument("-in2", "--infilename2", help="The file name of the second spatial map.")
    parser.add_argument("-nx", "--numberOfReportingCellsX", type=int, default=280, help="Number of reporting cells in x direction.")
    parser.add_argument("-ny", "--numberOfReportingCellsY", type=int, default=120, help="Number of reporting cells in y direction.")
    # Define shrinking factor. Since the algorithm is too memory consuming on the full
    # data, an integer shrinking factor >= 2 should be selected.
    parser.add_argument("-sf", "--shrinkingFactor", type=int, default=2, help="Shrinking factor >= 2 required to make algorithm applicable.")

    cmdArgs = vars(parser.parse_args())

    calculateEMD(cmdArgs["infilename1"], cmdArgs["infilename2"], cmdArgs["numberOfReportingCellsX"], cmdArgs["numberOfReportingCellsY"], cmdArgs["shrinkingFactor"])

if __name__ == "__main__":
    emd()
