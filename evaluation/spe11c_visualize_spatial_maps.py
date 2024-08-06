# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
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
#np.set_printoptions(threshold=sys.maxsize)

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

    delimiter = ','

    csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=5)
    csvData[:,2] = np.around(csvData[:,2], decimals=5)
    mask = (abs(csvData[:, 1] - 2475) < 1e-4)
    csvData = csvData[mask, :]
    ind = np.lexsort((csvData[:,0], csvData[:,2]))
    csvData = csvData[ind]
    for i in np.arange(0, nY):
        p[i, :] = csvData[i*nX:(i+1)*nX, 3] if len(csvData[0]) > 3 else 0
        s[i, :] = csvData[i*nX:(i+1)*nX, 4] if len(csvData[0]) > 4 else 0
        mCO2[i, :] = csvData[i*nX:(i+1)*nX, 5] if len(csvData[0]) > 5 else 0
        mH2O[i, :] = csvData[i*nX:(i+1)*nX, 6] if len(csvData[0]) > 6 else 0
        rhoG[i, :] = csvData[i*nX:(i+1)*nX, 7] if len(csvData[0]) > 7 else 0
        rhoL[i, :] = csvData[i*nX:(i+1)*nX, 8] if len(csvData[0]) > 8 else 0
        tmCO2[i, :] = csvData[i*nX:(i+1)*nX, 9] if len(csvData[0]) > 9 else 0
        temp[i, :] = csvData[i*nX:(i+1)*nX, 10] if len(csvData[0]) > 10 else 0

    p[p < 1e0] = float('nan')
    rhoG[rhoG < 1e-5] = float('nan')
    rhoL[rhoL < 1e-5] = float('nan')
    rhoG[s < 1e-3] = float('nan')
    rhoL[s > 1 - 1e-3] = float('nan')
    mH2O[s < 1e-3] = float('nan')
    mCO2[s > 1 - 1e-3] = float('nan')
    rhoG[s == float('nan')] = float('nan')
    rhoL[s == float('nan')] = float('nan')
    mH2O[s == float('nan')] = float('nan')
    mCO2[s == float('nan')] = float('nan')
    return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

def plotColorMesh(fig, x, y, z, idx, name, vmin, vmax, pRows, pCols):
    ax = fig.add_subplot(pRows, pCols, 1 + idx)
    im = ax.pcolormesh(x, y, z, shading='flat', cmap='viridis', vmin=vmin, vmax=vmax)
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
    fig.tight_layout()#rect=[0, 0.03, 1, 0.95])


def visualizeSpatialMaps():
    """Visualize spatial maps for Case C of the 11th SPE CSP"""

    font = {'size' : 14}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This script visualizes the spatial maps "
                    "on an evenly spaced grid as required by the description."
    )
    parser.add_argument("-t", "--time", default="5",
                        help="The time in years at which the spatial maps should be evaluated. "
                             "Assumes that the files are named 'spe11c_spatial_map_<X>y.csv', "
                             "where X is the given time. Defaults to '5'.")

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    cmdArgs = vars(parser.parse_args())
    time = cmdArgs["time"]
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    xSpace = np.arange(0.0, 8.4e3 + 25.0, 5.0e1)
    ySpace = np.arange(0.0, 3.6e3 + 5.0, 3.0e1)
    x, y = np.meshgrid(xSpace, ySpace)
    nX = xSpace.size - 1
    nY = ySpace.size - 1

    if len(groups) == 1:
        fig = plt.figure(figsize=(14, 8))
    else:
        figP = plt.figure(figsize=(14, 8))
        figS = plt.figure(figsize=(14, 8))
        figMCO2 = plt.figure(figsize=(14, 8))
        figMH2O = plt.figure(figsize=(14, 8))
        figRhoG = plt.figure(figsize=(14, 8))
        figRhoL = plt.figure(figsize=(14, 8))
        figTmCO2 = plt.figure(figsize=(14, 8))
        figTemp = plt.figure(figsize=(14, 8))

    if len(groups) == 1:
        pRows = 3
        pCols = 3
    elif len(groups) < 3:
        pRows = 1
        pCols = 2
    elif len(groups) < 5:
        pRows = 2
        pCols = 2
    elif len(groups) < 7:
        pRows = 2
        pCols = 3
    elif len(groups) < 10:
        pRows = 3
        pCols = 3
    elif len(groups) < 13:
        pRows = 3
        pCols = 4
    else:
        pRows = 4
        pCols = 4

    for i, group in zip(range(len(groups)), groups):
        if groupFolders:
            baseFolder = groupFolders[i]

        if group[-2] != '-':
            if not groupFolders:
                baseFolder = os.path.join(folder, group.lower())
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-2].lower(), f'result{group[-1]}')

        fileName = os.path.join(baseFolder, f'spe11c_spatial_map_{time}y.csv')
        p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileName, nX, nY)

        if len(groups) == 1:
            # scale pressure to bars
            plotColorMesh(fig, x, y, 1e-5*p, 0, "pressure [bar]", 200, 350, pRows, pCols)
            plotColorMesh(fig, x, y, s, 1, "gas saturation [-]", 0, 1, pRows, pCols)
            # scale mass fractions to g/kg
            plotColorMesh(fig, x, y, 1e3*mCO2, 2, "CO2 mass frac in liquid [g/kg]", 0, 70, pRows, pCols)
            plotColorMesh(fig, x, y, 1e3*mH2O, 3, "H2O mass frac in gas [g/kg]", 1, 4, pRows, pCols)
            plotColorMesh(fig, x, y, rhoG, 4, "gas phase density [kg/m3]", 0.8e3, 1.0e3, pRows, pCols)
            plotColorMesh(fig, x, y, rhoL, 5, "liquid phase density [kg/m3]", 0.99e3, 1.03e3, pRows, pCols)
            # scale mass to kilotons
            plotColorMesh(fig, x, y, 1e-6*tmCO2, 6, "total CO2 mass [kt]", 0, 4.5, pRows, pCols)
            plotColorMesh(fig, x, y, temp, 7, "temperature [°C]", 20, 70, pRows, pCols)
        else:
            # scale pressure to bars
            plotColorMesh(figP, x, y, 1e-5*p, i, group, 200, 350, pRows, pCols)
            plotColorMesh(figS, x, y, s, i, group, 0, 1, pRows, pCols)
            # scale mass fractions to g/kg
            plotColorMesh(figMCO2, x, y, 1e3*mCO2, i, group, 0, 70, pRows, pCols)
            plotColorMesh(figMH2O, x, y, 1e3*mH2O, i, group, 1, 4, pRows, pCols)
            plotColorMesh(figRhoG, x, y, rhoG, i, group, 0.8e3, 1.0e3, pRows, pCols)
            plotColorMesh(figRhoL, x, y, rhoL, i, group, 0.99e3, 1.03e3, pRows, pCols)
            # scale mass to kilotons
            plotColorMesh(figTmCO2, x, y, 1e-6*tmCO2, i, group, 0, 4.5, pRows, pCols)
            plotColorMesh(figTemp, x, y, temp, i, group, 20, 70, pRows, pCols)
    
    if len(groups) == 1:
        fig.suptitle(f'{groups[0]} at {time} years')
        fig.savefig(f'spe11c_{groups[0].lower()}_{time}y.png', bbox_inches='tight')
        print('File spe11c_' + f'{groups[0].lower()}_{time}y.png has been generated.')
    else:
        figP.suptitle(f'pressure [bar] at {time} years')
        figP.savefig(f'spe11c_pressure_{time}y.png', bbox_inches='tight')
        figS.suptitle(f'gas saturation [-] at {time} years')
        figS.savefig(f'spe11c_saturation_{time}y.png', bbox_inches='tight')
        figMCO2.suptitle(f'CO2 mass fraction in liquid [g/kg] at {time} years')
        figMCO2.savefig(f'spe11c_mco2_{time}y.png', bbox_inches='tight')
        figMH2O.suptitle(f'H2O mass fraction in gas [g/kg] at {time} years')
        figMH2O.savefig(f'spe11c_mh2o_{time}y.png', bbox_inches='tight')
        figRhoG.suptitle(f'gas phase density [kg/m3] at {time} years')
        figRhoG.savefig(f'spe11c_rhog_{time}y.png', bbox_inches='tight')
        figRhoL.suptitle(f'liquid phase density [kg/m3] at {time} years')
        figRhoL.savefig(f'spe11c_rhol_{time}y.png', bbox_inches='tight')
        figTmCO2.suptitle(f'total CO2 mass [kt] at {time} years')
        figTmCO2.savefig(f'spe11c_tmco2_{time}y.png', bbox_inches='tight')
        figTemp.suptitle(f'temperature [°C] at {time} years')
        figTemp.savefig(f'spe11c_temp_{time}y.png', bbox_inches='tight')
        print('Files spe11c_{pressure, saturation, mco2, mh2o, rhog, rhol, tmco2, temp}' + f'_{time}y.png have been generated.')

if __name__ == "__main__":
    visualizeSpatialMaps()
