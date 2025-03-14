# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3

""""
Script to reformat a result for the 11th SPE CSP
according to the reporting requirements
"""

import os
import argparse
import zipfile
import tempfile
import numpy as np
from is_notebook import is_notebook

def reformatSparseData(folder, case):
    fileName = os.path.join(folder, f"spe11{case}_time_series.csv")

    if os.path.exists(fileName):
        succeeded = True

        skip_header = 0
        with open(fileName, "r") as file:
            if not (file.readline()[0]).isnumeric():
                skip_header = 1

        csvData = np.genfromtxt(fileName, delimiter=",", skip_header=skip_header)

        if len(csvData.shape) == 1:
            print(f"  Couldn't detect any columns in {os.path.basename(fileName)}. Did you use the correct delimiter ','?")            
            succeeded = False
        else:
            columns = csvData.shape[1]

            if case == "a":
                requiredColumns = 13
            else:
                requiredColumns = 14

            if columns != requiredColumns:
                print(f"  The number of columns in {os.path.basename(fileName)} is {columns}, but should be {requiredColumns}.")
                if columns < requiredColumns:
                    print(f"  If a quantity can't be reported, the respective column should be filled with entries 'n/a' or 'nan'.")
                else:
                    print("  Please delete superfluous columns.")
                succeeded = False

            if csvData[0][0] != 0:
                if csvData[0][0] == 1.0:
                    csvData[0][0] = 0
                else:
                    print(f"  The first reported time step in {os.path.basename(fileName)} is {csvData[0][0]}, not the expected value 0.")
                    succeeded = False

            # four significant digits as required by the description
            fmt = ["%.3e"]*requiredColumns
            if case == "a":
                finalTime = 4.32e5
                header = "t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], mobB [kg], immB [kg], dissB [kg], sealB [kg], mC [m^2], sealTot [kg]"
                # change according to https://connect.spe.org/discussion/feedback-after-participant-workshop
                fmt[1] = fmt[2] = "%.5e"
            else:
                finalTime = 3.1536e10
                header = "t [s], p1 [Pa], p2 [Pa], mobA [kg], immA [kg], dissA [kg], sealA [kg], mobB [kg], immB [kg], dissB [kg], sealB [kg], mC [m^2], sealTot [kg], boundary [kg]"
                # change according to https://connect.spe.org/discussion/feedback-after-participant-workshop
                fmt[0] = "%.4e"

            if abs(csvData[-1][0] - finalTime)/finalTime > 1e-6:
                print(f"  The last reported time step in {os.path.basename(fileName)} is {csvData[-1][0]}, not the expected value {finalTime}.")
                succeeded = False

            for col in range(requiredColumns):
                if all(val < 0 for val in csvData[:, col]):
                    print(f'Only negative values in column {col}, replacing by nan.')
                    csvData[:, col] = np.nan

            newName = os.path.join(folder, "reformatted", f"spe11{case}_time_series.csv")
            np.savetxt(newName, csvData, fmt=fmt, delimiter=", ", header=header)
    else:
        print(f"  No sparse data file {os.path.basename(fileName)} exists.")
        print(f"  Did you use the correct filename and not introduce subfolders?")
        succeeded = False

    return succeeded


def reformatDenseData(folder, case):
    succeeded = True

    if case == "a":
        timeSteps = np.arange(0, 121, 1)
        unit = "h"
        numX = 280
        numY = 120
        numCells = numX*numY
        cellWidth = 1e-2
    elif case == "b":
        timeSteps = np.arange(0, 1001, 5)
        unit = "y"
        numX = 840
        numY = 120
        numCells = numX*numY
        cellWidth = 1e1
    else:
        timeSteps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
        unit = "y"
        numX = 168
        numY = 100
        numZ = 120
        numCells = numX*numY*numZ
        cellWidth = 50
        cellHeight = 10

    x = np.zeros(numCells)
    y = np.zeros(numCells)

    if case == "c":
        z = np.zeros(numCells)

        for k in range(numZ):
            for j in range(numY):
                x[j*numX:(j+1)*numX] = np.arange(0.5*cellWidth, cellWidth*numX, cellWidth)
                y[j*numX:(j+1)*numX] = cellWidth*(j + 0.5)
    else:
        for j in range(numY):
            x[j*numX:(j+1)*numX] = np.arange(0.5*cellWidth, cellWidth*numX, cellWidth)
            y[j*numX:(j+1)*numX] = cellWidth*(j + 0.5)

    missingSteps = []
    for time in timeSteps:
        fileName = os.path.join(folder, f"spe11{case}_spatial_map_{time}{unit}.csv")

        if not os.path.exists(fileName):
            missingSteps.append(time)
            succeeded = False

    if len(timeSteps) == len(missingSteps):
        print(f"  No dense data files spe11{case}_spatial_map_*{unit}.csv found in folder {folder}.")
        print(f"  Did you use the correct filenames and not introduce subfolders?")
    elif len(missingSteps) > 0:
        print(f"  No dense data files spe11{case}_spatial_map_*{unit}.csv found for times {missingSteps} {unit}.")

    for time in timeSteps:
        fileName = os.path.join(folder, f"spe11{case}_spatial_map_{time}{unit}.csv")
        if os.path.exists(fileName):
            skip_header = 0
            with open(fileName, "r") as file:
                if not (file.readline()[0]).isnumeric():
                    skip_header = 1

            csvData = np.genfromtxt(fileName, delimiter=',', skip_header=skip_header)

            if len(csvData.shape) == 1:
                print(f"  Couldn't detect any columns in {os.path.basename(fileName)}. Did you use the correct delimiter ','?")
                succeeded = False
            else:
                if case == "a" or case == "b":
                    csvData[:,0] = np.around(csvData[:,0], decimals=5)
                    csvData[:,1] = np.around(csvData[:,1], decimals=5)
                    ind = np.lexsort((csvData[:,0], csvData[:,1]))
                else:
                    csvData[:,0] = np.around(csvData[:,0], decimals=5)
                    csvData[:,1] = np.around(csvData[:,1], decimals=3)
                    csvData[:,2] = np.around(csvData[:,2], decimals=5)
                    ind = np.lexsort((csvData[:,0], csvData[:,1], csvData[:,2]))
                csvData = csvData[ind]
                if case == "a" or case == "b":
                    csvData[:,0] = x
                    csvData[:,1] = y

                requiredRows = numCells
                if case == "a":
                    requiredColumns = 9
                    header = "x [m], z [m], pressure [Pa], saturation [-], mCO2 [-], mH2O [-], rhoG [kg/m3], rhoL [kg/m3], tmCO2 [kg]"
                elif case == "b":
                    requiredColumns = 10
                    header = "x [m], z [m], pressure [Pa], saturation [-], mCO2 [-], mH2O [-], rhoG [kg/m3], rhoL [kg/m3], tmCO2 [kg], temp [°C]"
                else:
                    requiredColumns = 11
                    header = "x [m], y [m], z [m], pressure [Pa], saturation [-], mCO2 [-], mH2O [-], rhoG [kg/m3], rhoL [kg/m3], tmCO2 [kg], temp [°C]"

                rows = csvData.shape[0]
                columns = csvData.shape[1]

                if (rows != requiredRows):
                    print(f"  The number of reporting cells in {os.path.basename(fileName)} is {rows}, but should be {requiredRows}.")
                    succeeded = False

                if (columns != requiredColumns):
                    print(f"  The number of reported columns in {os.path.basename(fileName)} is {columns}, but should be {requiredColumns}.")
                    if columns < requiredColumns:
                        print(f"  If a quantity can't be reported, the respective column should be filled with entries 'n/a' or 'nan'.")
                        succeeded = False
                    else:
                        csvData = csvData[:, :requiredColumns]
                        print("  Superfluous columns deleted.")

                for col in range(requiredColumns):
                    if all(val < 0 for val in csvData[:, col]):
                        print(f'Only negative values in column {col}, replacing by nan.')
                        csvData[:, col] = np.nan

                newName = os.path.join(folder, "reformatted", f"spe11{case}_spatial_map_{time}{unit}.csv")
                np.savetxt(newName, csvData, fmt="%.3e", delimiter=", ", header=header)
        else:
            print(f"  Can't analyze in more detail, since initial file {os.path.basename(fileName)} is missing.")

    return succeeded


def reformat():
    """Reformat a result submitted to the 11th SPE CSP according to the requirements"""

    parser = argparse.ArgumentParser(
        description="This script reformats result submitted to the 11th SPE CSP according to the requirements."
    )
    parser.add_argument("-f", "--folder", default=".",
                        help="path to the folder containing the files for the result submission.")
    parser.add_argument("-c", "--case", default="A",
                        help="the CSP case, select one of {A, B, C}.")

    cmdArgs = vars(parser.parse_args())
    folder = cmdArgs["folder"]
    case = cmdArgs["case"].lower()

    print(f"\nReformatting for SPE11{case.upper()} in folder {folder}.\n")

    zipFileName = os.path.join(folder, f"spe11{case}.zip")

    createdTempDir = False
    if os.path.exists(zipFileName):
        tempDir = tempfile.TemporaryDirectory()
        createdTempDir = True
        folder = tempDir.name
        print(f"Found a zip file {zipFileName}. Unzipping to a temporary folder {folder}.\n")
        with zipfile.ZipFile(zipFileName, 'r') as zipRef:
            zipRef.extractall(folder)

    if not os.path.exists(os.path.join(folder, "reformatted")):
        os.makedirs(os.path.join(folder, "reformatted"))

    if reformatSparseData(folder, case):
        print("Successfully reformatted the sparse data.\n")
    else:
        print("The sparse data does not meet the requirements.\n")

    if reformatDenseData(folder, case):
        print("Successfully reformatted the dense data.\n")
    else:
        print("The dense data does not meet the requirements.\n")

    if createdTempDir:
        print(f"\nRemoving temporary folder {folder}.")
        tempDir.cleanup()

if __name__ == "__main__":
    reformat()
