#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from groups_and_colors import groups_and_colors

font = {'size' : 14}
matplotlib.rc('font', **font)

groups = ["CSIRO", "GEOS", "IFPEN", "OpenGoSim", "OPM", "Pau-INRIA", "PetroSim", "SLB", "UT-CSEE-PGE"]

numGroups = len(groups)
numSteps = 12

distances = np.loadtxt("distances.csv", delimiter=",")

fig, axs = plt.subplots(1, 1, figsize=(6, 3))

medianDistances = np.zeros((numSteps, numGroups), dtype=float)

for i in np.arange(0, numSteps):
    # The calculated distances have the unit of normalized mass times meter.
    # Multiply by 4.59, the injected mass of CO2 in g, and 100, to convert to g.cm.
    A = 459*distances[i*numGroups:(i+1)*numGroups, i*numGroups:(i+1)*numGroups]
    medianDistances[i, :] = np.median(A, axis=0) 

hours = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]
for group, i in zip(groups, np.arange(0, numGroups)):
    axs.plot(hours, medianDistances[:, i], color=groups_and_colors[group.lower()], label=group)

axs.set_xlabel(r'time [h]')
axs.set_ylabel(r'median EMD [g.cm]')
axs.set_ylim(0, 100)
#axs.set_title(r"CSP 11A Wasserstein distances over time")
axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig.savefig(f"spe11a_emd_over_time.png", bbox_inches='tight', dpi=300)
