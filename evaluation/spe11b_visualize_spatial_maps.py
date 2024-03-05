#!/usr/bin/env python3

""""
Script to visualize the spatial maps
on an evenly spaced grid as required by the description
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
np.set_printoptions(threshold=sys.maxsize)

def getFieldValues(fileName, nX, nY):
    p = np.zeros([nY, nX])
    s = np.zeros([nY, nX])
    mCO2 = np.zeros([nY, nX])
    mH2O = np.zeros([nY, nX])
    rhoG = np.zeros([nY, nX])
    rhoL = np.zeros([nY, nX])
    tmCO2 = np.zeros([nY, nX])
    temp = np.zeros([nY, nX])

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning 0 values.')
        return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

    skip_header = 0
    with open(fileName, "r") as file:
        if not (file.readline()[0]).isnumeric():
            skip_header = 1

    csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)
    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=5)
    ind = np.lexsort((csvData[:,0], csvData[:,1]))
    csvData = csvData[ind]
    for i in np.arange(0, nY):
        p[i, :] = csvData[i*nX:(i+1)*nX, 2] if len(csvData[0]) > 2 else 0
        s[i, :] = csvData[i*nX:(i+1)*nX, 3] if len(csvData[0]) > 3 else 0
        mCO2[i, :] = csvData[i*nX:(i+1)*nX, 4] if len(csvData[0]) > 4 else 0
        mH2O[i, :] = csvData[i*nX:(i+1)*nX, 5] if len(csvData[0]) > 5 else 0
        rhoG[i, :] = csvData[i*nX:(i+1)*nX, 6] if len(csvData[0]) > 6 else 0
        rhoL[i, :] = csvData[i*nX:(i+1)*nX, 7] if len(csvData[0]) > 7 else 0
        tmCO2[i, :] = csvData[i*nX:(i+1)*nX, 8] if len(csvData[0]) > 8 else 0
        temp[i, :] = csvData[i*nX:(i+1)*nX, 9] if len(csvData[0]) > 9 else 0

    return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

def plotColorMesh(fig, x, y, z, idx, name):
    ax = fig.add_subplot(331 + idx)
    im = ax.pcolormesh(x, y, z, shading='flat', cmap='coolwarm')
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    ax.axis('scaled')
    ax.set_title(f'{name}')
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbformat = matplotlib.ticker.ScalarFormatter()
    cbformat.set_powerlimits((-2,2))
    fig.colorbar(im, cax=cax, orientation='vertical', format=cbformat)


def visualizeSpatialMaps():
    """Visualize spatial maps for the SPE11 CSP"""

    font = {'size' : 12}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This script visualizes the spatial maps "
                    "on an evenly spaced grid as required by the description."
    )
    parser.add_argument("-t", "--time", default="5",
                        help="The time in years at which the spatial maps should be evaluated. "
                             "Assumes that the files are named 'spe11b_spatial_map_<X>y.csv', "
                             "where X is the given time. Defaults to '5'.")
    parser.add_argument('-g','--groups', nargs='+', help='<Required> names of groups', required=True)

    cmdArgs = vars(parser.parse_args())
    time = cmdArgs["time"]
    groups = cmdArgs["groups"]

    xSpace = np.arange(0.0, 8.4e3 + 5.0, 1.0e1)
    ySpace = np.arange(0.0, 3.6e3 + 5.0, 3.0e1)
    x, y = np.meshgrid(xSpace, ySpace)
    nX = xSpace.size - 1
    nY = ySpace.size - 1

    if len(groups) == 1:
        fig = plt.figure(figsize=(16, 7))
    else:
        figP = plt.figure(figsize=(16, 7))
        figS = plt.figure(figsize=(16, 7))
        figMCO2 = plt.figure(figsize=(16, 7))
        figMH2O = plt.figure(figsize=(16, 7))
        figRhoG = plt.figure(figsize=(16, 7))
        figRhoL = plt.figure(figsize=(16, 7))
        figTmCO2 = plt.figure(figsize=(16, 7))
        figTemp = plt.figure(figsize=(16, 7))

    for i, group in zip(range(len(groups)), groups):
        fileName = f'/media/bernd/bernd/spe11/{group.lower()}/spe11b_spatial_map_{time}y.csv'

        p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileName, nX, nY)
        if group.lower() == "ifpen": # did not report rhoL
            tmCO2 = rhoL
            rhoL = np.zeros([nY, nX])

        if len(groups) == 1:
            plotColorMesh(fig, x, y, p/1e5, 0, "pressure [bar]")
            plotColorMesh(fig, x, y, s, 1, "gas saturation [-]")
            plotColorMesh(fig, x, y, mCO2, 2, "CO2 mass frac in liquid [-]")
            plotColorMesh(fig, x, y, mH2O, 3, "H2O mass frac in gas [-]")
            plotColorMesh(fig, x, y, rhoG, 4, "gas phase density [kg/m3]")
            plotColorMesh(fig, x, y, rhoL, 5, "liquid phase density [kg/m3]")
            plotColorMesh(fig, x, y, tmCO2, 6, "total CO2 mass [kg]")
            plotColorMesh(fig, x, y, temp, 7, "temperature [°C]")
        else:
            plotColorMesh(figP, x, y, p/1e5, i, group)
            plotColorMesh(figS, x, y, s, i, group)
            plotColorMesh(figMCO2, x, y, mCO2, i, group)
            plotColorMesh(figMH2O, x, y, mH2O, i, group)
            plotColorMesh(figRhoG, x, y, rhoG, i, group)
            plotColorMesh(figRhoL, x, y, rhoL, i, group)
            plotColorMesh(figTmCO2, x, y, tmCO2, i, group)
            plotColorMesh(figTemp, x, y, temp, i, group)
    
    if len(groups) == 1:
        fig.suptitle(f'{groups[0]} at {time} years')
        fig.savefig(f'spe11b_{groups[0].lower()}_{time}y.png', bbox_inches='tight')
        print('File spe11b_' + f'{groups[0].lower()}_{time}y.png has been generated.')
    else:
        figP.suptitle(f'pressure [bar] at {time} years')
        figP.savefig(f'spe11b_pressure_{time}y.png', bbox_inches='tight')
        figS.suptitle(f'gas saturation [-] at {time} years')
        figS.savefig(f'spe11b_saturation_{time}y.png', bbox_inches='tight')
        figMCO2.suptitle(f'CO2 mass fraction in liquid [-] at {time} years')
        figMCO2.savefig(f'spe11b_mco2_{time}y.png', bbox_inches='tight')
        figMH2O.suptitle(f'H2O mass fraction in gas [-] at {time} years')
        figMH2O.savefig(f'spe11b_mh2o_{time}y.png', bbox_inches='tight')
        figRhoG.suptitle(f'gas phase density [kg/m3] at {time} years')
        figRhoG.savefig(f'spe11b_rhog_{time}y.png', bbox_inches='tight')
        figRhoL.suptitle(f'liquid phase density [kg/m3] at {time} years')
        figRhoL.savefig(f'spe11b_rhol_{time}y.png', bbox_inches='tight')
        figTmCO2.suptitle(f'total CO2 mass [kg] at {time} years')
        figTmCO2.savefig(f'spe11b_tmco2_{time}y.png', bbox_inches='tight')
        figTemp.suptitle(f'temperature [°C] at {time} years')
        figTemp.savefig(f'spe11b_temp_{time}y.png', bbox_inches='tight')
        print('Files spe11b_{pressure, saturation, mco2, mh2o, rhog, rhol, tmco2, temp}' + f'_{time}y.png have been generated.')

if __name__ == "__main__":
    visualizeSpatialMaps()
