# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import numpy as np
import os.path
import argparse
import emd

def calculateWassersteinDistances():
    """Calculate Wasserstein distances for Case A of the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This calculates Wasserstein distances based on given spatial maps."
    )

    parser.add_argument('-g','--groups', nargs='+', help='names of groups', required=True)

    parser.add_argument('-gf','--groupfolders', nargs='+', help='paths to group folders', required=False)

    parser.add_argument('-f','--folder', help='path to folder containing group subfolders', required=False)

    parser.add_argument('--hours', nargs='+', help='hour(s) for which to calculate the distances', type=int, required=False, default=24)

    cmdArgs = vars(parser.parse_args())
    groups = cmdArgs["groups"]
    groupFolders = cmdArgs["groupfolders"]
    folder = cmdArgs["folder"]
    hours = cmdArgs["hours"]

    # number of reporting cells in each coordinate direction
    nX = 280
    nY = 120
    shrinkingFactor = 2

    numGroups = len(groups)
    numHours = len(hours)
    distances = np.zeros((numGroups*numHours, numGroups*numHours))

    for iHour, hourI in zip(range(numHours), hours):
        for jHour, hourJ in zip(range(numHours), hours):
            if iHour > jHour: continue

            for iGroup, groupI in zip(range(numGroups), groups):
                if groupFolders:
                    baseFolderI = groupFolders[iGroup]
                else:
                    if groupI[-2] != '-':
                        baseFolderI = os.path.join(folder, groupI.lower())
                    else:
                        baseFolderI = os.path.join(folder, groupI[:-2].lower(), f'result{groupI[-1]}')

                fileNameI = os.path.join(baseFolderI, f'spe11a_spatial_map_{hourI}h.csv')
                if (os.path.exists(fileNameI)):
                    row = iHour*numGroups + iGroup

                    for jGroup, groupJ in zip(range(numGroups), groups):
                        if jGroup <= iGroup and hourJ == hourI: continue

                        if groupFolders:
                            baseFolderJ = groupFolders[jGroup]
                        else:
                            if groupJ[-2] != '-':
                                baseFolderJ = os.path.join(folder, groupJ.lower())
                            else:
                                baseFolderJ = os.path.join(folder, groupJ[:-2].lower(), f'result{groupJ[-1]}')

                        fileNameJ = os.path.join(baseFolderJ, f'spe11a_spatial_map_{hourJ}h.csv')
                        if (os.path.exists(fileNameJ)):
                            col = jHour*numGroups + jGroup

                            distances[row][col] = emd.calculateEMD(fileNameI, fileNameJ, nX, nY, shrinkingFactor)

                            print(f'{hourI}, {hourJ}, {groupI}, {groupJ} -> ({row}, {col}): {distances[row][col]}')

    distances = distances + distances.T - np.diag(distances.diagonal())

    np.savetxt("distances.csv", distances, fmt='%.3e', delimiter=",")
    print("Distances have been stored in distances.csv.")
    print("The calculated distances have the unit of normalized mass times meter.")

if __name__ == "__main__":
    calculateWassersteinDistances()
