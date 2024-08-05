import numpy as np
import os.path
import emd


baseFileNames = ['/media/bernd/bernd/spe11/early/grayscale/spe11a_csiro_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_geos_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_ifpen_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_opengosim_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_opm_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_pau-inria_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_petrosim_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_slb_grayscale_',
                 '/media/bernd/bernd/spe11/early/grayscale/spe11a_ut-csee-pge_grayscale_']

numGroups = len(baseFileNames)
distances = np.zeros((numGroups*12, numGroups*12))

for hourI in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]:
    for hourJ in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]:
        if hourJ != hourI: continue

        for i, baseFileNameI in zip(range(numGroups), baseFileNames):
            fileNameI = baseFileNameI + str(hourI) + 'h.png'
            if (os.path.exists(fileNameI)):
                row = int((hourI/10 - 1)*numGroups + i)

                for j, baseFileNameJ in zip(range(numGroups), baseFileNames):
                    if j <= i and hourJ == hourI: continue

                    fileNameJ = baseFileNameJ + str(hourJ) + 'h.png'
                    if (os.path.exists(fileNameJ)):
                        col = int((hourJ/10 - 1)*numGroups + j)

                        distances[row][col] = emd.calculateEMD(fileNameI, fileNameJ)

                        print(f'{hourI}, {hourJ}, {i}, {j} -> ({row}, {col}): {distances[row][col]}')

distances = distances + distances.T - np.diag(distances.diagonal())

np.savetxt("distances.csv", distances, delimiter=",")
