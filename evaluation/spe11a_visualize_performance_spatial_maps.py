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
from is_notebook import is_notebook

def getFieldValues(fileName, nX, nY):
    cvol = np.zeros([nY, nX])
    arat = np.zeros([nY, nX])
    co2_max_norm_res = np.zeros([nY, nX])
    h2o_max_norm_res = np.zeros([nY, nX])
    co2_mb_error = np.zeros([nY, nX])
    h2o_mb_error = np.zeros([nY, nX])
    post_est = np.zeros([nY, nX])

    if os.path.isfile(fileName):
        print(f'Processing {fileName}.')
    else:
        print(f'No file {fileName} found. Returning 0 values.')
        return cvol, arat, co2_max_norm_res, h2o_max_norm_res, co2_mb_error, h2o_mb_error, post_est

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
        cvol[i, :] = csvData[i*nX:(i+1)*nX, 2] if len(csvData[0]) > 2 else 0
        arat[i, :] = csvData[i*nX:(i+1)*nX, 3] if len(csvData[0]) > 3 else 0
        co2_max_norm_res[i, :] = csvData[i*nX:(i+1)*nX, 4] if len(csvData[0]) > 4 else 0
        h2o_max_norm_res[i, :] = csvData[i*nX:(i+1)*nX, 5] if len(csvData[0]) > 5 else 0
        co2_mb_error[i, :] = csvData[i*nX:(i+1)*nX, 6] if len(csvData[0]) > 6 else 0
        h2o_mb_error[i, :] = csvData[i*nX:(i+1)*nX, 7] if len(csvData[0]) > 7 else 0
        post_est[i, :] = csvData[i*nX:(i+1)*nX, 8] if len(csvData[0]) > 8 else 0

    return cvol, arat, co2_max_norm_res, h2o_max_norm_res, co2_mb_error, h2o_mb_error, post_est

def plotColorMesh(fig, x, y, z, idx, name, pRows, pCols):
    ax = fig.add_subplot(pRows, pCols, 1 + idx)
    im = ax.pcolormesh(x, y, z, shading='flat', cmap='viridis')
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


def visualizePerformanceSpatialMaps():
    """Visualize performance spatial maps for Case A of the 11th SPE CSP"""

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
    groups = [x.lower() for x in cmdArgs["groups"]]
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
        figCVol = plt.figure(figsize=(14, 8))
        figARat = plt.figure(figsize=(14, 8))
        figCO2MNR = plt.figure(figsize=(14, 8))
        figH2OMNR = plt.figure(figsize=(14, 8))
        figCO2MBE = plt.figure(figsize=(14, 8))
        figH2OMBE = plt.figure(figsize=(14, 8))
        figPostEst = plt.figure(figsize=(14, 8))

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

        if not group[-1].isnumeric():
            if not groupFolders:
                baseFolder = os.path.join(folder, group, 'spe11a')
        else:
            if not groupFolders:
                baseFolder = os.path.join(folder, group[:-1], 'spe11a', f'result{group[-1]}')

        fileName = os.path.join(baseFolder, f'spe11a_performance_spatial_map_{time}h.csv')
        cvol, arat, co2_max_norm_res, h2o_max_norm_res, co2_mb_error, h2o_mb_error, post_est = getFieldValues(fileName, nX, nY)

        if len(groups) == 1:
            plotColorMesh(fig, x, y, 1e6*cvol, 0, r"cell volume [cm$^3$]", pRows, pCols)
            plotColorMesh(fig, x, y, arat, 1, "aspect ratio [-]", pRows, pCols)
            plotColorMesh(fig, x, y, co2_max_norm_res, 2, "CO2 maximum normalized residual [-]", pRows, pCols)
            plotColorMesh(fig, x, y, h2o_max_norm_res, 3, "H2O maximum normalized residual [-]", pRows, pCols)
            plotColorMesh(fig, x, y, co2_mb_error, 4, "CO2 mass-balance error [-]", pRows, pCols)
            plotColorMesh(fig, x, y, h2o_mb_error, 5, "H2O mass-balance error [-]", pRows, pCols)
            plotColorMesh(fig, x, y, post_est, 6, "a posteriori error estimate [-]", pRows, pCols)
        else:
            plotColorMesh(figCVol, x, y, 1e6*cvol, i, group, pRows, pCols)
            plotColorMesh(figARat, x, y, arat, i, group, pRows, pCols)
            plotColorMesh(figCO2MNR, x, y, co2_max_norm_res, i, group, pRows, pCols)
            plotColorMesh(figH2OMNR, x, y, h2o_max_norm_res, i, group, pRows, pCols)
            plotColorMesh(figCO2MBE, x, y, co2_mb_error, i, group, pRows, pCols)
            plotColorMesh(figH2OMBE, x, y, h2o_mb_error, i, group, pRows, pCols)
            plotColorMesh(figPostEst, x, y, post_est, i, group, pRows, pCols)
    
    if len(groups) == 1:
        fig.suptitle(f'{groups[0]} at {time} hours')
        if not is_notebook():
            fig.savefig(f'spe11a_{groups[0].lower()}_performance+_{time}h.png', bbox_inches='tight')
            print('File spe11a_' + f'{groups[0].lower()}_performance_{time}h.png has been generated.')
    else:
        figCVol.suptitle(r'cell volume [cm$^3$] at' + f' {time} hours')
        figARat.suptitle(f'aspect ratio [-] at {time} hours')
        figCO2MNR.suptitle(f'CO2 maximum normalized residual [-] at {time} hours')
        figH2OMNR.suptitle(f'H2O maximum normalized residual [-] at {time} hours')
        figCO2MBE.suptitle(f'CO2 mass-balance error [-] at {time} hours')
        figH2OMBE.suptitle(f'H2O mass-balance error [-] at {time} hours')
        figPostEst.suptitle(f'a posteriori error estimate [-] at {time} hours')
        if not is_notebook():
            figCVol.savefig(f'spe11a_cvol_{time}h.png', bbox_inches='tight')
            figARat.savefig(f'spe11a_arat_{time}h.png', bbox_inches='tight')
            figCO2MNR.savefig(f'spe11a_co2_max_norm_res_{time}h.png', bbox_inches='tight')
            figH2OMNR.savefig(f'spe11a_h2o_max_norm_res_{time}h.png', bbox_inches='tight')
            figCO2MBE.savefig(f'spe11a_co2_mb_error_{time}h.png', bbox_inches='tight')
            figH2OMBE.savefig(f'spe11a_h2o_mb_error_{time}h.png', bbox_inches='tight')
            figPostEst.savefig(f'spe11a_post_est_{time}h.png', bbox_inches='tight')
            print('Files spe11a_{cvol, arat, {co2, h2o}_max_norm_res, {co2, h2o}_mb_error, post_est}' + f'_{time}h.png have been generated.')

if __name__ == "__main__":
    visualizePerformanceSpatialMaps()
