#!/usr/bin/env python3

""""
Script to generate a grayscale image of the combined gas saturation
and CO2 concentration values on an evenly spaced grid as required
by the benchmark description
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def getFieldValues(fileName, nX, nY):
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
    if "geos" in fileName:
        skip_header = 2

    csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)
    if "geos" in fileName: # remove additional y coordinate column
        csvData = np.delete(csvData, 1, 1) 
    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=5)
    ind = np.lexsort((csvData[:,0], csvData[:,1]))
    csvData = csvData[ind]
    idx = 8 if not "ifpen" in fileName else 7
    for i in np.arange(0, nY):
        tmCO2[i, :] = csvData[i*nX:(i+1)*nX, idx] if len(csvData[0]) > idx else 0

    tmCO2[np.isnan(tmCO2)] = 0

    return tmCO2

def plotColorMesh(fig, x, y, z, outFileName):
    ax = fig.gca()
    vmin = 0.0
    vmax = 1.8

    im = ax.pcolormesh(x, y, z, shading='flat', cmap='gist_gray')#, vmin=vmin, vmax=vmax)
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    ax.axis('scaled')
    plt.axis('off')
    
    fig.savefig(outFileName, bbox_inches='tight', dpi=96, pad_inches=0)
    print(f'File {outFileName} has been generated.')


def generateGrayScale():
    """Generate a grayscale a spatial map for the FluidFlower benchmark"""

    groups = ["csiro", "geos", "ifpen", "opengosim", "opm", "pau-inria", "petrosim", "slb", "ut-csee-pge"]

    xSpace = np.arange(0.0, 2.8 + 5.0e-3, 1.0e-2)
    ySpace = np.arange(0.0, 1.2 + 5.0e-3, 1.0e-2)
    x, y = np.meshgrid(xSpace, ySpace)
    nX = xSpace.size - 1
    nY = ySpace.size - 1

    # The target size of the resulting image is 280x120 pixels, one pixel per reporting cell.
    # The added pixels correspond to the automatically added whitespace which can be removed
    # only when saving the figure to disk.
    my_dpi = 96
    fig = plt.figure(figsize=((280+82)/my_dpi, (120+37)/my_dpi), dpi=my_dpi)

    for group in groups:
        for hour in np.arange(0, 121):
            inFileName = f"/media/bernd/bernd/spe11/early/{group}/spe11a_spatial_map_{hour}h.csv"
            if "ut-csee" in group:
                inFileName = f"/media/bernd/bernd/spe11/early/{group}/result1/spe11a_spatial_map_{hour}h.csv"

            tmCO2 = getFieldValues(inFileName, nX, nY)

            outFileName = f"/media/bernd/bernd/spe11/early/grayscale/spe11a_{group}_grayscale_{hour}h.png"

            plotColorMesh(fig, x, y, tmCO2, outFileName)

if __name__ == "__main__":
    generateGrayScale()
