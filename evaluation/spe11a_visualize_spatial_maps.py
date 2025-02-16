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
import seaborn as sns
import groups_and_colors
from is_notebook import is_notebook
from round_to_digits import round_to_digits

def getFieldValues(fileName, nX, nY):
    p = np.zeros([nY, nX]); p[:] = np.nan
    s = np.zeros([nY, nX]); s[:] = np.nan
    mCO2 = np.zeros([nY, nX]); mCO2[:] = np.nan
    mH2O = np.zeros([nY, nX]); mH2O[:] = np.nan
    rhoG = np.zeros([nY, nX]); rhoG[:] = np.nan
    rhoL = np.zeros([nY, nX]); rhoL[:] = np.nan
    tmCO2 = np.zeros([nY, nX]); tmCO2[:] = np.nan

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning nans.')
        return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2

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
        # round to four significant digits as required by the description
        p[i, :] = round_to_digits(csvData[i*nX:(i+1)*nX, 2], 4) if len(csvData[0]) > 2 else 0
        s[i, :] = round_to_digits(csvData[i*nX:(i+1)*nX, 3], 4) if len(csvData[0]) > 3 else 0
        mCO2[i, :] = round_to_digits(csvData[i*nX:(i+1)*nX, 4], 4) if len(csvData[0]) > 4 else 0
        mH2O[i, :] = round_to_digits(csvData[i*nX:(i+1)*nX, 5], 4) if len(csvData[0]) > 5 else 0
        rhoG[i, :] = round_to_digits(csvData[i*nX:(i+1)*nX, 6], 4) if len(csvData[0]) > 6 else 0
        rhoL[i, :] = round_to_digits(csvData[i*nX:(i+1)*nX, 7], 4) if len(csvData[0]) > 7 else 0
        tmCO2[i, :] = round_to_digits(csvData[i*nX:(i+1)*nX, 8], 4) if len(csvData[0]) > 8 else 0

    p[p < 1e0] = np.nan
    rhoG[rhoG < 1e-5] = np.nan
    rhoL[rhoL < 1e-5] = np.nan
    rhoG[s < 1e-3] = -1
    rhoL[s > 1 - 1e-3] = np.nan
    mH2O[s < 1e-3] = -1
    mCO2[s > 1 - 1e-3] = np.nan
    rhoG[np.isnan(s)] = np.nan
    rhoL[np.isnan(s)] = np.nan
    mH2O[np.isnan(s)] = np.nan
    mCO2[np.isnan(s)] = np.nan
    return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2

def plotColorMesh(fig, x, y, z, idx, name, pRows, pCols, cmap='viridis', vmin=None, vmax=None):
    if cmap == 'viridis' or cmap == 'coolwarm':
        resetUnder = False
    else:
        resetUnder = True

    if isinstance(cmap, str):
        cmap = matplotlib.colormaps[cmap]
        
    if resetUnder:
        cmap.set_under([1, 1, 1])
    cmap.set_bad([0.5, 0.5, 0.5])

    if vmin is None:
        vmin = np.nanmin(np.where(z > 0, z, np.inf))

    if vmax is None:
        vmax = np.nanmax(z)

    ax = fig.add_subplot(pRows, pCols, 1 + idx)
    if vmax == vmin:
        im = ax.pcolormesh(x, y, z, shading='flat', cmap=cmap)
    else:
        im = ax.pcolormesh(x, y, z, shading='flat', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    ax.axis('scaled')
    ax.set_title(f'{name}')
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbformat = matplotlib.ticker.ScalarFormatter()
    cbformat.set_powerlimits((-2,2))
#    fig.colorbar(im, cax=cax, orientation='vertical', format=cbformat)
#    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.subplots_adjust(right=1.0)
    cbar_ax = fig.add_axes([1.02, 0.14, 0.025, 0.7])
    fig.colorbar(im, cax=cbar_ax, format=cbformat)


def visualizeSpatialMaps():
    """Visualize spatial maps for Case A of the 11th SPE CSP"""

    font = {'size' : 14}
    matplotlib.rc('font', **font)

    parser = argparse.ArgumentParser(
        description="This script visualizes the spatial maps "
                    "on an evenly spaced grid as required by the description."
    )
    parser.add_argument("-t", "--time", default="24",
                        help="The time in hours at which the spatial maps should be evaluated. "
                             "Assumes that the files are named 'spe11a_spatial_map_<X>h.csv', "
                             "where X is the given time. Defaults to '24'.")

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    cmdArgs = vars(parser.parse_args())
    time = cmdArgs["time"]
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]

    xSpace = np.arange(0.0, 2.8 + 5.0e-3, 1.0e-2)
    ySpace = np.arange(0.0, 1.2 + 5.0e-3, 1.0e-2)
    x, y = np.meshgrid(xSpace, ySpace)
    nX = xSpace.size - 1
    nY = ySpace.size - 1

    if len(groups) == 1:
        fig = plt.figure(figsize=(14, 8))
    else:
        figP = plt.figure(figsize=(15, 9))
        figS = plt.figure(figsize=(15, 9))
        figMCO2 = plt.figure(figsize=(15, 9))
        figMH2O = plt.figure(figsize=(15, 9))
        figRhoG = plt.figure(figsize=(15, 9))
        figRhoL = plt.figure(figsize=(15, 9))
        figTmCO2 = plt.figure(figsize=(14, 7.2))

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

    # select file that contains impermeable cells with 'nan' pressure values
    fileNameSLB = os.path.join(folder, 'slb', 'spe11a', 'result1', f'spe11a_spatial_map_{time}h.csv')
    pSLB, s, mCO2, mH2O, rhoG, rhoL, tmCO2 = getFieldValues(fileNameSLB, nX, nY)

    for i, group in zip(range(len(groups)), groups):
        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group.lower(), 'spe11a')
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1].lower(), 'spe11a', f'result{group[-1]}')

        fileName = os.path.join(baseFolder, f'spe11a_spatial_map_{time}h.csv')
        p, s, mCO2, mH2O, rhoG, rhoL, tmCO2 = getFieldValues(fileName, nX, nY)
        p[np.isnan(pSLB)] = np.nan
        s[np.isnan(pSLB)] = np.nan
        mCO2[np.isnan(pSLB)] = np.nan
        mH2O[np.isnan(pSLB)] = np.nan
        rhoG[np.isnan(pSLB)] = np.nan
        rhoL[np.isnan(pSLB)] = np.nan
        tmCO2[np.isnan(pSLB)] = np.nan

        cmap = groups_and_colors.mass_cmap

        if len(groups) == 1:
            # scale pressure to bars
            plotColorMesh(fig, x, y, 1e-5*p, 0, "pressure [bar]", pRows, pCols)
            plotColorMesh(fig, x, y, s, 1, "gas saturation [-]", pRows, pCols, cmap)
            # scale mass fractions to g/kg
            plotColorMesh(fig, x, y, 1e3*mCO2, 2, "CO2 mass frac in liquid [g/kg]", pRows, pCols, cmap)
            plotColorMesh(fig, x, y, 1e3*mH2O, 3, "H2O mass frac in gas [g/kg]", pRows, pCols, 'icefire')
            plotColorMesh(fig, x, y, rhoG, 4, "gas phase density [kg/m3]", pRows, pCols, 'icefire')
            plotColorMesh(fig, x, y, rhoL, 5, "liquid phase density [kg/m3]", pRows, pCols, 'icefire')
            # scale mass to grams
            plotColorMesh(fig, x, y, 1e3*tmCO2, 6, "total CO2 mass [g]", pRows, pCols, cmap)
        else:
            # scale pressure to bars
            plotColorMesh(figP, x, y, 1e-5*p, i, group, pRows, pCols, 'viridis', 1.1, 1.22)
            plotColorMesh(figS, x, y, s, i, group, pRows, pCols, cmap, 0, 1)
            # scale mass fractions to g/kg
            plotColorMesh(figMCO2, x, y, 1e3*mCO2, i, group, pRows, pCols, cmap, 0, 2)
            plotColorMesh(figMH2O, x, y, 1e3*mH2O, i, group, pRows, pCols, 'icefire', 8.1, 8.7)
            plotColorMesh(figRhoG, x, y, rhoG, i, group, pRows, pCols, 'icefire', 2.0, 2.2)
            plotColorMesh(figRhoL, x, y, rhoL, i, group, pRows, pCols, 'icefire', 9.973e2, 9.987e2)
            # scale mass to grams
            plotColorMesh(figTmCO2, x, y, 1e3*tmCO2, i, group, pRows, pCols, cmap, 0, 1e-3)
    
    if len(groups) == 1:
        fig.suptitle(f'{groups[0]} at {time} hours')
        if not is_notebook():
            fig.savefig(f'spe11a_{groups[0].lower()}_{time}h.png', bbox_inches='tight')
            print('File spe11a_' + f'{groups[0].lower()}_{time}h.png has been generated.')
    else:
        figP.suptitle(f'pressure [bar] at {time} hours')
        figS.suptitle(f'gas saturation [-] at {time} hours')
        figMCO2.suptitle(f'CO2 mass fraction in liquid [g/kg] at {time} hours')
        figMH2O.suptitle(f'H2O mass fraction in gas [g/kg] at {time} hours')
        figRhoG.suptitle(f'gas phase density [kg/m3] at {time} hours')
        figRhoL.suptitle(f'liquid phase density [kg/m3] at {time} hours')
        figTmCO2.suptitle(f'                          total CO2 mass [g] at {time} hours')

        figP.savefig(f'spe11a_pressure_{time}h.png', bbox_inches='tight')
        figS.savefig(f'spe11a_saturation_{time}h.png', bbox_inches='tight')
        figMCO2.savefig(f'spe11a_mco2_{time}h.png', bbox_inches='tight')
        figMH2O.savefig(f'spe11a_mh2o_{time}h.png', bbox_inches='tight')
        figRhoG.savefig(f'spe11a_rhog_{time}h.png', bbox_inches='tight')
        figRhoL.savefig(f'spe11a_rhol_{time}h.png', bbox_inches='tight')
        figTmCO2.savefig(f'spe11a_tmco2_{time}h.png', bbox_inches='tight', dpi=300)
        print('Files spe11a_{pressure, saturation, mco2, mh2o, rhog, rhol, tmco2}' + f'_{time}h.png have been generated.')

if __name__ == "__main__":
    visualizeSpatialMaps()
