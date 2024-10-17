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

def getFieldValues(fileName, nX, nY, nZ):
    p = np.zeros([nX, nY, nZ]); p[:] = np.nan
    s = np.zeros([nX, nY, nZ]); s[:] = np.nan
    mCO2 = np.zeros([nX, nY, nZ]); mCO2[:] = np.nan
    mH2O = np.zeros([nX, nY, nZ]); mH2O[:] = np.nan
    rhoG = np.zeros([nX, nY, nZ]); rhoG[:] = np.nan
    rhoL = np.zeros([nX, nY, nZ]); rhoL[:] = np.nan
    tmCO2 = np.zeros([nX, nY, nZ]); tmCO2[:] = np.nan
    temp = np.zeros([nX, nY, nZ]); temp[:] = np.nan

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning nans.')
        return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

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
        print(f'Cannot find unique x ({len(np.unique(csvData[:,0]))} instead of {nX}) or y ({len(np.unique(csvData[:,1]))} instead of {nY}) coordinates. Returning nans.')
        return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

    for i in np.arange(0, nX):
        for j in np.arange(0, nY):
            p[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 3] if len(csvData[0]) > 3 else float('nan')
            s[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 4] if len(csvData[0]) > 4 else float('nan')
            mCO2[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 5] if len(csvData[0]) > 5 else float('nan')
            mH2O[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 6] if len(csvData[0]) > 6 else float('nan')
            rhoG[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 7] if len(csvData[0]) > 7 else float('nan')
            rhoL[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 8] if len(csvData[0]) > 8 else float('nan')
            tmCO2[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 9] if len(csvData[0]) > 9 else float('nan')
            temp[i, j, :] = csvData[(i*nY + j)*nZ:(i*nY + j + 1)*nZ, 10] if len(csvData[0]) > 10 else float('nan')

    p[p < 1e0] = float('nan')
    rhoG[rhoG < 1e-5] = float('nan')
    rhoL[rhoL < 1e-5] = float('nan')
    rhoG[s < 1e-3] = float('nan')
    rhoL[s > 1 - 1e-3] = float('nan')
    mH2O[s < 1e-3] = float('nan')
    mCO2[s > 1 - 1e-3] = float('nan')
    rhoG[np.isnan(s)] = float('nan')
    rhoL[np.isnan(s)] = float('nan')
    mH2O[np.isnan(s)] = float('nan')
    mCO2[np.isnan(s)] = float('nan')
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

    parser.add_argument('-cp','--cutplane', help='the plane where to cut, out of {uv, uw, vw}', required=False, default='uw')

    parser.add_argument('-idx','--cutindex', help='index where to cut', required=False, type=int, default=49)

    cmdArgs = vars(parser.parse_args())
    time = cmdArgs["time"]
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    cutPlane = cmdArgs["cutplane"]
    cutIndex = cmdArgs["cutindex"]

    if len(groups) == 1:
        fig = plt.figure(figsize=(14, 8))
    else:
        figP = plt.figure(figsize=(15, 9))
        figS = plt.figure(figsize=(15, 9))
        figMCO2 = plt.figure(figsize=(15, 9))
        figMH2O = plt.figure(figsize=(15, 9))
        figRhoG = plt.figure(figsize=(15, 9))
        figRhoL = plt.figure(figsize=(15, 9))
        figTmCO2 = plt.figure(figsize=(15, 9))
        figTemp = plt.figure(figsize=(15, 9))

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
    elif len(groups) < 17:
        pRows = 4
        pCols = 4
    elif len(groups) < 21:
        pRows = 4
        pCols = 5
    else:
        pRows = 5
        pCols = 5

    nX = 168
    nY = 100
    nZ = 120

    # select file that contains impermeable cells with 'nan' pressure values
    fileNameTetra = os.path.join(folder, 'tetratech-rps', 'spe11c', 'result2', f'spe11c_spatial_map_{time}y.csv')
    p, s, mCO2Tetra, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileNameTetra, nX, nY, nZ)

    for i, group in zip(range(len(groups)), groups):
        if groupFolders:
            baseFolder = groupFolders[i]

        if group[-2] != '-':
            if not groupFolders:
                baseFolder = os.path.join(folder, group.lower(), 'spe11c')
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-2].lower(), 'spe11c', f'result{group[-1]}')

        fileName = os.path.join(baseFolder, f'spe11c_spatial_map_{time}y.csv')
        p3, s3, mCO23, mH2O3, rhoG3, rhoL3, tmCO23, temp3 = getFieldValues(fileName, nX, nY, nZ)
        p3[np.isnan(mCO2Tetra)] = float('nan')
        s3[np.isnan(mCO2Tetra)] = float('nan')
        mCO23[np.isnan(mCO2Tetra)] = float('nan')
        mH2O3[np.isnan(mCO2Tetra)] = float('nan')
        rhoG3[np.isnan(mCO2Tetra)] = float('nan')
        rhoL3[np.isnan(mCO2Tetra)] = float('nan')
        tmCO23[np.isnan(mCO2Tetra)] = float('nan')

        if cutPlane == 'uw':
            xSpace = np.arange(0.0, 8.4e3 + 25.0, 5.0e1)
            ySpace = np.arange(0.0, 3.6e3 + 15.0, 3.0e1)
            p = np.transpose(p3[0:nX, cutIndex, 0:nZ])
            s = np.transpose(s3[0:nX, cutIndex, 0:nZ])
            mCO2 = np.transpose(mCO23[0:nX, cutIndex, 0:nZ])
            mH2O = np.transpose(mH2O3[0:nX, cutIndex, 0:nZ])
            rhoG = np.transpose(rhoG3[0:nX, cutIndex, 0:nZ])
            rhoL = np.transpose(rhoL3[0:nX, cutIndex, 0:nZ])
            tmCO2 = np.transpose(tmCO23[0:nX, cutIndex, 0:nZ])
            temp = np.transpose(temp3[0:nX, cutIndex, 0:nZ])
        elif cutPlane == 'vw':
            xSpace = np.arange(0.0, 5.0e3 + 25.0, 5.0e1)
            ySpace = np.arange(0.0, 1.2e3 + 5.0, 1.0e1)
            p = np.transpose(p3[cutIndex, 0:nY, 0:nZ])
            s = np.transpose(s3[cutIndex, 0:nY, 0:nZ])
            mCO2 = np.transpose(mCO23[cutIndex, 0:nY, 0:nZ])
            mH2O = np.transpose(mH2O3[cutIndex, 0:nY, 0:nZ])
            rhoG = np.transpose(rhoG3[cutIndex, 0:nY, 0:nZ])
            rhoL = np.transpose(rhoL3[cutIndex, 0:nY, 0:nZ])
            tmCO2 = np.transpose(tmCO23[cutIndex, 0:nY, 0:nZ])
            temp = np.transpose(temp3[cutIndex, 0:nY, 0:nZ])

        x, y = np.meshgrid(xSpace, ySpace)
        if cutPlane == 'vw':
            for j, v in zip(range(0, len(x[0])), x[0]):
                y[:, j] = 2*(y[:, j] + 150*(1 - np.square((v - 2500)/2500)) + v/500)

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
        figP.savefig(f'spe11c_{cutPlane}_{cutIndex}_pressure_{time}y.png', bbox_inches='tight')
        figS.suptitle(f'gas saturation [-] at {time} years')
        figS.savefig(f'spe11c_{cutPlane}_{cutIndex}_saturation_{time}y.png', bbox_inches='tight')
        figMCO2.suptitle(f'CO2 mass fraction in liquid [g/kg] at {time} years')
        figMCO2.savefig(f'spe11c_{cutPlane}_{cutIndex}_mco2_{time}y.png', bbox_inches='tight')
        figMH2O.suptitle(f'H2O mass fraction in gas [g/kg] at {time} years')
        figMH2O.savefig(f'spe11c_{cutPlane}_{cutIndex}_mh2o_{time}y.png', bbox_inches='tight')
        figRhoG.suptitle(f'gas phase density [kg/m3] at {time} years')
        figRhoG.savefig(f'spe11c_{cutPlane}_{cutIndex}_rhog_{time}y.png', bbox_inches='tight')
        figRhoL.suptitle(f'liquid phase density [kg/m3] at {time} years')
        figRhoL.savefig(f'spe11c_{cutPlane}_{cutIndex}_rhol_{time}y.png', bbox_inches='tight')
        figTmCO2.suptitle(f'total CO2 mass [kt] at {time} years')
        figTmCO2.savefig(f'spe11c_{cutPlane}_{cutIndex}_tmco2_{time}y.png', bbox_inches='tight')
        figTemp.suptitle(f'temperature [°C] at {time} years')
        figTemp.savefig(f'spe11c_{cutPlane}_{cutIndex}_temp_{time}y.png', bbox_inches='tight')
        print(f'Files spe11c_{cutPlane}_{cutIndex}_' + '{pressure, saturation, mco2, mh2o, rhog, rhol, tmco2, temp}' + f'_{time}y.png have been generated.')

if __name__ == "__main__":
    visualizeSpatialMaps()
