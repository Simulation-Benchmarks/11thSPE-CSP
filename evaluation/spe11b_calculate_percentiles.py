# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from groups_and_colors import groups_and_colors
from is_notebook import is_notebook

def assembleTimeSeries():
    """Calculate percentiles for Case B of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script calculates the percentiles of the time series quantities "
                    "required by the CSP description."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups, taking reported numbers', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('-c','--calculated', nargs='+', help='names of groups, taking calculated numbers for Boxes A and B')

    parser.add_argument('-t','--tablefolder', help='path to folder containing calculated tables')

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    calculated = []
    if cmdArgs["calculated"]:
        calculated = cmdArgs["calculated"]
        groups = sorted(groups + calculated)
        tableFolder = cmdArgs["tablefolder"]

    if calculated:
        mobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_mobA_from_spatial_maps.csv')
        mobileFromSpatialMapsA = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
        tSpatialMaps = mobileFromSpatialMapsA['time_s']/60/60
        # project initial values to 1e-1 hours for improved visualization
        tSpatialMaps[0] = 1e-1
        immobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_immA_from_spatial_maps.csv')
        immobileFromSpatialMapsA = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
        dissolvedFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_dissA_from_spatial_maps.csv')
        dissolvedFromSpatialMapsA = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
        sealFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_sealA_from_spatial_maps.csv')
        sealFromSpatialMapsA = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

        mobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_mobB_from_spatial_maps.csv')
        mobileFromSpatialMapsB = np.genfromtxt(mobileFromSpatialMapsFileName, delimiter=',', names=True)
        immobileFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_immB_from_spatial_maps.csv')
        immobileFromSpatialMapsB = np.genfromtxt(immobileFromSpatialMapsFileName, delimiter=',', names=True)
        dissolvedFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_dissB_from_spatial_maps.csv')
        dissolvedFromSpatialMapsB = np.genfromtxt(dissolvedFromSpatialMapsFileName, delimiter=',', names=True)
        sealFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_sealB_from_spatial_maps.csv')
        sealFromSpatialMapsB = np.genfromtxt(sealFromSpatialMapsFileName, delimiter=',', names=True)

        convectionFromSpatialMapsFileName = os.path.join(tableFolder, 'spe11b_mC_from_spatial_maps.csv')
        convectionFromSpatialMaps = np.genfromtxt(convectionFromSpatialMapsFileName, delimiter=',', names=True)

    numGroups = len(groups)
    numTimeSteps = 10001
    t = np.linspace(0, 1000*365*24*60*60, num=numTimeSteps)
    p1 = np.zeros((numTimeSteps, numGroups));
    p2 = np.zeros((numTimeSteps, numGroups));
    mobA = np.zeros((numTimeSteps, numGroups));
    immA = np.zeros((numTimeSteps, numGroups));
    dissA = np.zeros((numTimeSteps, numGroups));
    sealA = np.zeros((numTimeSteps, numGroups));
    mobB = np.zeros((numTimeSteps, numGroups));
    immB = np.zeros((numTimeSteps, numGroups));
    dissB = np.zeros((numTimeSteps, numGroups));
    sealB = np.zeros((numTimeSteps, numGroups));
    mC = np.zeros((numTimeSteps, numGroups));
    sealTot = np.zeros((numTimeSteps, numGroups));

    for i, group in zip(range(numGroups), groups):
        if groupFolders:
            baseFolder = groupFolders[i]
        else:
            if not group[-1].isnumeric():
                baseFolder = os.path.join(folder, group.lower(), 'spe11b')
            else:
                baseFolder = os.path.join(folder, group[:-1].lower(), 'spe11b', f'result{group[-1]}')

        fileName = os.path.join(baseFolder, 'spe11b_time_series.csv')
        print(f'Processing {fileName}.')

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        delimiter = ','

        csvData = np.genfromtxt(fileName, delimiter=delimiter, skip_header=skip_header)

        if len(csvData[:, 0]) != numTimeSteps:
            print(f'Interpolating from {len(csvData[:, 0])} time steps.')
            p1[:, i] = np.interp(t, csvData[:, 0], csvData[:, 1])
            p2[:, i] = np.interp(t, csvData[:, 0], csvData[:, 2])
            sealTot[:, i] = np.interp(t, csvData[:, 0], csvData[:, 12])
        else:
            p1[:, i] = csvData[:, 1]
            p2[:, i] = csvData[:, 2]
            sealTot[:, i] = csvData[:, 12]

        if group in calculated:
            print(f'Interpolating from {mobileFromSpatialMapsA['time_s']} time steps.')
            columnName = group.lower().replace('-', '')
            mobA[:, i] = np.interp(t, mobileFromSpatialMapsA['time_s'], mobileFromSpatialMapsA[columnName])
            immA[:, i] = np.interp(t, immobileFromSpatialMapsA['time_s'], immobileFromSpatialMapsA[columnName])
            dissA[:, i] = np.interp(t, dissolvedFromSpatialMapsA['time_s'], dissolvedFromSpatialMapsA[columnName])
            sealA[:, i] = np.interp(t, sealFromSpatialMapsA['time_s'], sealFromSpatialMapsA[columnName])
            mobB[:, i] = np.interp(t, mobileFromSpatialMapsB['time_s'], mobileFromSpatialMapsB[columnName])
            immB[:, i] = np.interp(t, immobileFromSpatialMapsB['time_s'], immobileFromSpatialMapsB[columnName])
            dissB[:, i] = np.interp(t, dissolvedFromSpatialMapsB['time_s'], dissolvedFromSpatialMapsB[columnName])
            sealB[:, i] = np.interp(t, sealFromSpatialMapsB['time_s'], sealFromSpatialMapsB[columnName])
            mC[:, i] = np.interp(t, convectionFromSpatialMaps['time_s'], convectionFromSpatialMaps[columnName])
        elif len(csvData[:, 0]) != numTimeSteps:
            mobA[:, i] = np.interp(t, csvData[:, 0], csvData[:, 3])
            immA[:, i] = np.interp(t, csvData[:, 0], csvData[:, 4])
            dissA[:, i] = np.interp(t, csvData[:, 0], csvData[:, 5])
            sealA[:, i] = np.interp(t, csvData[:, 0], csvData[:, 6])
            mobB[:, i] = np.interp(t, csvData[:, 0], csvData[:, 7])
            immB[:, i] = np.interp(t, csvData[:, 0], csvData[:, 8])
            dissB[:, i] = np.interp(t, csvData[:, 0], csvData[:, 9])
            sealB[:, i] = np.interp(t, csvData[:, 0], csvData[:, 10])
            mC[:, i] = np.interp(t, csvData[:, 0], csvData[:, 11])
        else:
            mobA[:, i] = csvData[:, 3]
            immA[:, i] = csvData[:, 4]
            dissA[:, i] = csvData[:, 5]
            sealA[:, i] = csvData[:, 6]
            mobB[:, i] = csvData[:, 7]
            immB[:, i] = csvData[:, 8]
            dissB[:, i] = csvData[:, 9]
            sealB[:, i] = csvData[:, 10]
            mC[:, i] = csvData[:, 11]


    numPercentiles = 21
    p1Percentiles = np.zeros((numTimeSteps, numPercentiles+1))
    p2Percentiles = np.zeros((numTimeSteps, numPercentiles+1))
    mobAPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    immAPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    dissAPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    sealAPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    mobBPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    immBPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    dissBPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    sealBPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    mCPercentiles = np.zeros((numTimeSteps, numPercentiles+1))
    sealTotPercentiles = np.zeros((numTimeSteps, numPercentiles+1))

    pHeader = 'time [s]'
    massHeader = 'time [s]'
    mCHeader = 'time [s]'
    p1Percentiles[:, 0] = t
    p2Percentiles[:, 0] = t
    mobAPercentiles[:, 0] = t
    immAPercentiles[:, 0] = t
    dissAPercentiles[:, 0] = t
    sealAPercentiles[:, 0] = t
    mobBPercentiles[:, 0] = t
    immBPercentiles[:, 0] = t
    dissBPercentiles[:, 0] = t
    sealBPercentiles[:, 0] = t
    mCPercentiles[:, 0] = t
    sealTotPercentiles[:, 0] = t
    for i, p in zip(range(numPercentiles), np.arange(0, 101, 5)):
        p1Percentiles[:, i+1] = np.nanpercentile(p1, p, axis=1)
        p2Percentiles[:, i+1] = np.nanpercentile(p2, p, axis=1)
        mobAPercentiles[:, i+1] = np.nanpercentile(mobA, p, axis=1)
        immAPercentiles[:, i+1] = np.nanpercentile(immA, p, axis=1)
        dissAPercentiles[:, i+1] = np.nanpercentile(dissA, p, axis=1)
        sealAPercentiles[:, i+1] = np.nanpercentile(sealA, p, axis=1)
        mobBPercentiles[:, i+1] = np.nanpercentile(mobB, p, axis=1)
        immBPercentiles[:, i+1] = np.nanpercentile(immB, p, axis=1)
        dissBPercentiles[:, i+1] = np.nanpercentile(dissB, p, axis=1)
        sealBPercentiles[:, i+1] = np.nanpercentile(sealB, p, axis=1)
        mCPercentiles[:, i+1] = np.nanpercentile(mC, p, axis=1)
        sealTotPercentiles[:, i+1] = np.nanpercentile(sealTot, p, axis=1)
        pHeader = pHeader + f', P{p} [Pa]' 
        massHeader = massHeader + f', P{p} [kg]' 
        mCHeader = mCHeader + f', P{p} [m]' 

    if not is_notebook():
        np.savetxt('spe11b_p1_percentiles.csv', p1Percentiles, fmt='%.5e', delimiter=', ', header=pHeader)
        np.savetxt('spe11b_p2_percentiles.csv', p2Percentiles, fmt='%.5e', delimiter=', ', header=pHeader)
        np.savetxt('spe11b_mobA_percentiles.csv', mobAPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_immA_percentiles.csv', immAPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_dissA_percentiles.csv', dissAPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_sealA_percentiles.csv', sealAPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_mobB_percentiles.csv', mobBPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_immB_percentiles.csv', immBPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_dissB_percentiles.csv', dissBPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_sealB_percentiles.csv', sealBPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)
        np.savetxt('spe11b_mC_percentiles.csv', mCPercentiles, fmt='%.5e', delimiter=', ', header=mCHeader)
        np.savetxt('spe11b_sealTot_percentiles.csv', sealTotPercentiles, fmt='%.5e', delimiter=', ', header=massHeader)

if __name__ == "__main__":
    assembleTimeSeries()
