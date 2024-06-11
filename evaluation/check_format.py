#!/usr/bin/env python3

""""
Script to check the format requirements of a
result to be submitted to the 11th SPE CSP
"""

import os
import argparse
import zipfile
import tempfile
import numpy as np

def checkSparseData(folder, case):
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
                print(f"  The first reported time step in {os.path.basename(fileName)} is {csvData[0][0]}, not the expected value 0.")
                succeeded = False

            if case == "a":
                finalTime = 4.32e5
            else:
                finalTime = 3.1536e10

            if abs(csvData[-1][0] - finalTime)/finalTime > 1e-6:
                print(f"  The last reported time step in {os.path.basename(fileName)} is {csvData[-1][0]}, not the expected value {finalTime}.")
                succeeded = False
    else:
        print(f"  No sparse data file {os.path.basename(fileName)} exists.")
        print(f"  Did you use the correct filename and not introduce subfolders?")
        succeeded = False

    return succeeded


def checkDenseData(folder, case):
    succeeded = True

    if case == "a":
        timeSteps = np.arange(0, 121, 1)
        unit = "h"
    elif case == "b":
        timeSteps = np.arange(0, 1001, 5)
        unit = "y"
    else:
        timeSteps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]
        unit = "y"

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

    fileName = os.path.join(folder, f"spe11{case}_spatial_map_0{unit}.csv")
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
            csvData[:,0] = np.around(csvData[:,0], decimals=5)
            csvData[:,1] = np.around(csvData[:,1], decimals=5)
            ind = np.lexsort((csvData[:,0], csvData[:,1]))
            csvData = csvData[ind]

            if case == "a":
                requiredRows = 280*120
                requiredColumns = 9
            elif case == "b":
                requiredRows = 840*120
                requiredColumns = 10
            else:
                requiredRows = 168*100*120
                requiredColumns = 11

            rows = csvData.shape[0]
            columns = csvData.shape[1]

            if (rows != requiredRows):
                print(f"  The number of reporting cells in {os.path.basename(fileName)} is {rows}, but should be {requiredRows}.")
                succeeded = False

            if (columns != requiredColumns):
                print(f"  The number of reported columns in {os.path.basename(fileName)} is {columns}, but should be {requiredColumns}.")
                succeeded = False
                if columns < requiredColumns:
                    print(f"  If a quantity can't be reported, the respective column should be filled with entries 'n/a' or 'nan'.")
                else:
                    print("  Please delete superfluous columns.")
    else:
        print(f"  Can't analyze in more detail, since initial file {os.path.basename(fileName)} is missing.")

    return succeeded


def checkQuestionnaire(folder, case):
    fileName = os.path.join(folder, f"spe11{case}_questionnaire.xlsx")

    if os.path.exists(fileName):
        return True
    else:
        print(f"  No questionnaire {os.path.basename(fileName)} exists.")
        return False

def checkFormat():
    """check the format requirements of a result to be submitted to the 11th SPE CSP"""

    parser = argparse.ArgumentParser(
        description="This script checks the format requirements of a result "
                    "to be submitted to the 11th SPE CSP."
    )
    parser.add_argument("-f", "--folder", default=".",
                        help="path to the folder containing the files for the result submission.")
    parser.add_argument("-c", "--case", default="A",
                        help="the CSP case, select one of {A, B, C}.")

    cmdArgs = vars(parser.parse_args())
    folder = cmdArgs["folder"]
    case = cmdArgs["case"].lower()

    print(f"\nChecking result for CSP 11{case.upper()} in folder {folder}.\n")

    zipFileName = os.path.join(folder, f"spe11{case}.zip")

    createdTempDir = False
    if os.path.exists(zipFileName):
        tempDir = tempfile.TemporaryDirectory()
        createdTempDir = True
        folder = tempDir.name
        print(f"Found a zip file {zipFileName}. Unzipping to a temporary folder {folder}.\n")
        with zipfile.ZipFile(zipFileName, 'r') as zipRef:
            zipRef.extractall(folder)

    if checkSparseData(folder, case):
        print("Successfully checked the sparse data.\n")
    else:
        print("The sparse data does not meet the requirements.\n")

    if checkDenseData(folder, case):
        print("Successfully checked the dense data.\n")
    else:
        print("The dense data does not meet the requirements.\n")

    if checkQuestionnaire(folder, case):
        print("Successfully checked the questionnaire.")
    else:
        print("The questionnaire does not meet the requirements.")

    if createdTempDir:
        print(f"\nRemoving temporary folder {folder}.")
        tempDir.cleanup()

if __name__ == "__main__":
    checkFormat()
