# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
""""
Script to visualize the spatial maps
on an evenly spaced grid as required by the description
"""

import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import groups_and_colors
from round_to_digits import round_to_digits
import utils

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

    try:
        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)
    except:
        print(f'Reading into table from file {fileName} failed.')
        return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

    csvData[:,0] = np.around(csvData[:,0], decimals=5)
    csvData[:,1] = np.around(csvData[:,1], decimals=3)
    csvData[:,2] = np.around(csvData[:,2], decimals=5)
    ind = np.lexsort((csvData[:,2], csvData[:,1], csvData[:,0]))
    csvData = csvData[ind]

    if len(np.unique(csvData[:,0])) != nX or len(np.unique(csvData[:,1])) != nY:
        print(f'Cannot find unique x ({len(np.unique(csvData[:,0]))} instead of {nX}) or y ({len(np.unique(csvData[:,1]))} instead of {nY}) coordinates. Returning nans.')
        return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

    p = np.reshape(round_to_digits(csvData[:, 3], 4), (nX, nY, nZ)) if len(csvData[0]) > 3 else np.nan
    s = np.reshape(round_to_digits(csvData[:, 4], 4), (nX, nY, nZ)) if len(csvData[0]) > 4 else np.nan
    mCO2 = np.reshape(round_to_digits(csvData[:, 5], 4), (nX, nY, nZ)) if len(csvData[0]) > 5 else np.nan
    mH2O = np.reshape(round_to_digits(csvData[:, 6], 4), (nX, nY, nZ)) if len(csvData[0]) > 6 else np.nan
    rhoG = np.reshape(round_to_digits(csvData[:, 7], 4), (nX, nY, nZ)) if len(csvData[0]) > 7 else np.nan
    rhoL = np.reshape(round_to_digits(csvData[:, 8], 4), (nX, nY, nZ)) if len(csvData[0]) > 8 else np.nan
    tmCO2 = np.reshape(round_to_digits(csvData[:, 9], 4), (nX, nY, nZ)) if len(csvData[0]) > 9 else np.nan
    temp = np.reshape(round_to_digits(csvData[:, 10], 4), (nX, nY, nZ)) if len(csvData[0]) > 10 else np.nan

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
    return p, s, mCO2, mH2O, rhoG, rhoL, tmCO2, temp

def plotColorMesh(fig, x, y, z, idx, name, pRows, pCols, cmap='viridis', vmin=None, vmax=None, title=None):
    if cmap == 'viridis' or cmap == 'coolwarm':
        resetUnder = False
    else:
        resetUnder = True

    if isinstance(cmap, str):
        cmap = matplotlib.colormaps[cmap]
        
    if resetUnder:
        cmap.set_under([1, 1, 1])
    cmap.set_bad([0.5, 0.5, 0.5])

    onecbar = True
    if vmin is None:
        vmin = np.nanmin(np.where(z > 0, z, np.inf))
        onecbar = False
    if vmax is None:
        vmax = np.nanmax(z)
        onecbar = False

    ax = fig.add_subplot(pRows, pCols, 1 + idx)
    if vmax == vmin:
        im = ax.pcolormesh(x, y, z, shading='flat', cmap=cmap)
        onecbar = False
    else:
        im = ax.pcolormesh(x, y, z, shading='flat', cmap=cmap, vmin=vmin, vmax=vmax)

    ax.axis([x.min(), x.max(), y.min(), y.max()])
    ax.axis('scaled')
    ax.set_title(f'{name}')
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.plot(x[0], y[-1], color='k')
    ax.plot(x[0], y[0], color='k')
#    ax.tick_params(axis="x", direction="in")
#    ax.set_yticks([600, 800])
#    ax.set_yticklabels([300, 400])
#    ax.tick_params(axis="y", direction="in", pad=-30)
#    ax.grid(visible=True, which='both', axis='y')
#    ax.plot(x[0], y[70], color='k', linewidth=0.5)

    if onecbar:
        fig.subplots_adjust(right=1.0)
        cbar_ax = fig.add_axes([1.03, 0.14, 0.025, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.ax.tick_params(axis='y', labelsize='large')
        cbar.ax.ticklabel_format(style='sci', axis='y', scilimits=(-1,2))
        if title:
            cbar.ax.set_ylabel(title, fontsize='large')
            cbar.ax.yaxis.set_label_position('left')
    else:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbformat = matplotlib.ticker.ScalarFormatter()
        cbformat.set_powerlimits((-2,2))
        fig.colorbar(im, cax=cax, orientation='vertical', format=cbformat)
        fig.tight_layout()


def visualizeSpatialMaps():
    """Visualize spatial maps for Case C of the 11th SPE CSP"""

    utils.set_fonts(size=14)

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
    nX = 168
    nY = 100
    nZ = 120

    if len(groups) == 1:
        pRows = 3
        pCols = 3
        figsize = (14, 7.5)
    elif len(groups) < 3:
        pRows = 1
        pCols = 2
        figsize = (14, 6.3)
    elif len(groups) < 5:
        pRows = 2
        pCols = 2
        if cutPlane == 'vw':
            figsize = (11, 7.2)
        else:
            figsize = (14, 7.2)
    elif len(groups) < 7:
        pRows = 2
        pCols = 3
        figsize = (14, 6.3)
    elif len(groups) < 10:
        pRows = 3
        pCols = 3
        if cutPlane == 'vw':
            figsize = (14, 9.5)
        else:
            figsize = (14, 8)
    elif len(groups) < 13:
        pRows = 3
        pCols = 4
        if cutPlane == 'vw':
            figsize = (14, 7.5)
        else:
            figsize = (14, 6.5)
    elif len(groups) < 17:
        pRows = 4
        pCols = 4
        if cutPlane == 'vw':
            figsize = (14, 10.5)
        else:
            figsize = (14, 8.5)
    elif len(groups) < 21:
        pRows = 4
        pCols = 5
        if cutPlane == 'vw':
            figsize = (14, 8)
        else:
            figsize = (14, 7.2)
    else:
        pRows = 5
        pCols = 5
        if cutPlane == 'vw':
            figsize = (14, 10)
        else:
            figsize = (14, 9)

    if len(groups) == 1:
        fig = plt.figure(figsize=figsize)
    else:
        figP = plt.figure(figsize=figsize)
        figS = plt.figure(figsize=figsize)
        figMCO2 = plt.figure(figsize=figsize)
        figMH2O = plt.figure(figsize=figsize)
        figRhoG = plt.figure(figsize=figsize)
        figRhoL = plt.figure(figsize=figsize)
        figTmCO2 = plt.figure(figsize=figsize)
        figTemp = plt.figure(figsize=figsize)

        titleP = f'pressure [Pa] at {time} years'
        titleS = f'gas saturation [-] at {time} years'
        titleMCO2 = f'CO$_2$ mass fraction in liquid [-] at {time} years'
        titleMH2O = f'H2O mass fraction in gas [-] at {time} years'
        titleRhoG = f'gas phase density [kg/m3] at {time} years'
        titleRhoL = f'liquid phase density [kg/m3] at {time} years'
        titleTmCO2 = f'total CO$_2$ mass [kg] at {time} years'
        titleTemp = f'temperature [°C] at {time} years'

    # select file that contains impermeable cells with 'nan' pressure values
    fileNameTetra = os.path.join(folder, 'tetratech', 'spe11c', 'result2', f'spe11c_spatial_map_{time}y.csv')
    p, s, mCO2Tetra, mH2O, rhoG, rhoL, tmCO2, temp = getFieldValues(fileNameTetra, nX, nY, nZ)

    for i, group in zip(range(len(groups)), groups):
        if groupFolders:
            baseFolder = groupFolders[i]

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group.lower(), 'spe11c')
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1].lower(), 'spe11c', f'result{group[-1]}')

        fileName = os.path.join(baseFolder, f'spe11c_spatial_map_{time}y.csv')
        p3, s3, mCO23, mH2O3, rhoG3, rhoL3, tmCO23, temp3 = getFieldValues(fileName, nX, nY, nZ)
        p3[np.isnan(mCO2Tetra)] = np.nan
        s3[np.isnan(mCO2Tetra)] = np.nan
        mCO23[np.isnan(mCO2Tetra)] = np.nan
        mH2O3[np.isnan(mCO2Tetra)] = np.nan
        rhoG3[np.isnan(mCO2Tetra)] = np.nan
        rhoL3[np.isnan(mCO2Tetra)] = np.nan
        tmCO23[np.isnan(mCO2Tetra)] = np.nan

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

        cmap = groups_and_colors.mass_cmap

        if len(groups) == 1:
            plotColorMesh(fig, x, y, p, 0, "pressure [Pa]", pRows, pCols)
            plotColorMesh(fig, x, y, s, 1, "gas saturation [-]", pRows, pCols, cmap)
            plotColorMesh(fig, x, y, mCO2, 2, "CO$_2$ mass frac in liquid [-]", pRows, pCols, cmap)
            plotColorMesh(fig, x, y, mH2O, 3, "H2O mass frac in gas [-]", pRows, pCols, 'icefire')
            plotColorMesh(fig, x, y, rhoG, 4, "gas phase density [kg/m3]", pRows, pCols, 'icefire')
            plotColorMesh(fig, x, y, rhoL, 5, "liquid phase density [kg/m3]", pRows, pCols, 'icefire')
            plotColorMesh(fig, x, y, tmCO2, 6, "total CO$_2$ mass [kg]", pRows, pCols, cmap)
            plotColorMesh(fig, x, y, temp, 7, "temperature [°C]", pRows, pCols, 'coolwarm')
        else:
            plotColorMesh(figP, x, y, p, i, group, pRows, pCols, 'viridis', 2e7, 3.5e7, titleP)
            plotColorMesh(figS, x, y, s, i, group, pRows, pCols, cmap, 0, 1, titleS)
            plotColorMesh(figMCO2, x, y, mCO2, i, group, pRows, pCols, cmap, 0, 7e-2, titleMCO2)
            plotColorMesh(figMH2O, x, y, mH2O, i, group, pRows, pCols, 'icefire', 1e-3, 4e-3, titleMH2O)
            plotColorMesh(figRhoG, x, y, rhoG, i, group, pRows, pCols, 'icefire', 0.8e3, 1.0e3, titleRhoG)
            plotColorMesh(figRhoL, x, y, rhoL, i, group, pRows, pCols, 'icefire', 0.99e3, 1.03e3, titleRhoL)
            plotColorMesh(figTmCO2, x, y, tmCO2, i, group, pRows, pCols, cmap, 0, 3e6, titleTmCO2)
            plotColorMesh(figTemp, x, y, temp, i, group, pRows, pCols, 'coolwarm', 20, 70, titleTemp)
    
    if len(groups) == 1:
        fig.suptitle(f'{groups[0]} at {time} years')
        fig.savefig(f'spe11c_{groups[0].lower()}_{time}y.png', bbox_inches='tight', dpi=300)
        print('File spe11c_' + f'{groups[0].lower()}_{time}y.png has been generated.')
    else:
        figP.subplots_adjust(wspace=0.03)
        figS.subplots_adjust(wspace=0.03)
        figMCO2.subplots_adjust(wspace=0.03)
        figMH2O.subplots_adjust(wspace=0.03)
        figRhoG.subplots_adjust(wspace=0.03)
        figRhoL.subplots_adjust(wspace=0.03)
        figTmCO2.subplots_adjust(wspace=0.03)
        figTemp.subplots_adjust(wspace=0.03)

        figP.savefig(f'spe11c_{cutPlane}_{cutIndex}_pressure_{time}y.png', bbox_inches='tight', dpi=300)
        figS.savefig(f'spe11c_{cutPlane}_{cutIndex}_saturation_{time}y.png', bbox_inches='tight', dpi=300)
        figMCO2.savefig(f'spe11c_{cutPlane}_{cutIndex}_mco2_{time}y.png', bbox_inches='tight', dpi=300)
        figMH2O.savefig(f'spe11c_{cutPlane}_{cutIndex}_mh2o_{time}y.png', bbox_inches='tight', dpi=300)
        figRhoG.savefig(f'spe11c_{cutPlane}_{cutIndex}_rhog_{time}y.png', bbox_inches='tight', dpi=300)
        figRhoL.savefig(f'spe11c_{cutPlane}_{cutIndex}_rhol_{time}y.png', bbox_inches='tight', dpi=300)
        figTmCO2.savefig(f'spe11c_{cutPlane}_{cutIndex}_tmco2_{time}y.png', bbox_inches='tight', dpi=300)
        figTemp.savefig(f'spe11c_{cutPlane}_{cutIndex}_temp_{time}y.png', bbox_inches='tight', dpi=300)
        print(f'Files spe11c_{cutPlane}_{cutIndex}_' + '{pressure, saturation, mco2, mh2o, rhog, rhol, tmco2, temp}' + f'_{time}y.png have been generated.')

if __name__ == "__main__":
    visualizeSpatialMaps()
