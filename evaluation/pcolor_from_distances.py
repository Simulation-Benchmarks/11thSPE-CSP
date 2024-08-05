import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

font = {'family' : 'normal',
        'weight' : 'normal',
        'size' : 14}
matplotlib.rc('font', **font)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "monospace",
})

groups = ["CSIRO", "GEOS", "IFPEN", "OpenGoSim", "OPM", "Pau-INRIA", "PetroSim", "SLB", "UT-CSEE-PGE"]

numGroups = len(groups)

distances = np.loadtxt("distances.csv", delimiter=",")

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","red"])

fig, axs = plt.subplots(1, 2, figsize=(12, 3))

# The calculated distances have the unit of normalized mass times meter.
# Multiply by 4.59, the injected mass of CO2 in g, and 100, to convert to g.cm.
A = 459*distances[:numGroups, :numGroups]
# remove values related to LANL
#A = np.delete(A, 5, 0)
#A = np.delete(A, 5, 1)
meanA = np.median(A, axis=0) 
A = A + np.diag(meanA)
A = np.flip(A, 0)

#groups.remove("LANL")
#colors.remove("#FF1493")
axs[0].pcolor(A, cmap=cmap)
axs[0].set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
axs[0].set_yticklabels(reversed(groups))
#for i in range(numGroups-1):
#    axs[0].get_yticklabels()[i].set_color(colors[numGroups-2-i])
axs[0].set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
axs[0].set_xticklabels(groups, rotation=45, ha="right")
#for i in range(numGroups-1):
#    axs[0].get_xticklabels()[i].set_color(colors[i])
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
})
for (i,j), value in np.ndenumerate(A):
    if i >= 8-j:
        if value > 0:
            if i == 8-j:
                axs[0].text(numGroups-i-0.5, numGroups-j-0.6, f'{int(value):03d}', style='italic', ha='center', va='center')
            else:
                axs[0].text(numGroups-i-0.5, numGroups-j-0.6, f'{int(value):03d}', ha='center', va='center')
        else:
            axs[0].text(numGroups-i-0.5, numGroups-j-0.6, '-', ha='center', va='center')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "monospace",
})
axs[0].set_title(r"\textrm{\textbf{10 hours}}")

A = 459*distances[11*numGroups:, 11*numGroups:]
# remove values related to LANL
#A = np.delete(A, 5, 0)
#A = np.delete(A, 5, 1)
meanA = np.median(A, axis=0) 
A = A + np.diag(meanA)
A = np.flip(A, 0)

axs[1].pcolor(A, cmap=cmap)
axs[1].tick_params(axis='y', which='both', left=False, right=True, labelleft=False, labelright=True)
axs[1].set_yticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
axs[1].set_yticklabels(reversed(groups))
#for i in range(numGroups-1):
#    axs[1].get_yticklabels()[i].set_color(colors[numGroups-2-i])
axs[1].set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
axs[1].set_xticklabels(groups, rotation=45, ha="right")
#for i in range(numGroups-1):
#    axs[1].get_xticklabels()[i].set_color(colors[i])
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
})
for (i,j), value in np.ndenumerate(A):
    if i >= 8-j:
        if value > 0:
            if i == 8-j:
                axs[1].text(numGroups-i-0.5, numGroups-j-0.6, f'{int(value):03d}', style='italic', ha='center', va='center')
            else:
                axs[1].text(numGroups-i-0.5, numGroups-j-0.6, f'{int(value):03d}', ha='center', va='center')
        else:
            axs[1].text(numGroups-i-0.5, numGroups-j-0.6, '-', ha='center', va='center')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "monospace",
})
axs[1].set_title(r"\textrm{\textbf{120 hours}}")

fig.savefig(f"pcolor_distances.png", bbox_inches='tight', dpi=300)
