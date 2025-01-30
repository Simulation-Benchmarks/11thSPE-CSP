# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import seaborn

groups_and_colors = {
    "calgary": "#90B0DD",
    "cau-kiel": "#C40F10",
    "csiro": "#0B5197",
    "ctc-cne": "#FFA052",
    "darts": "#FF5A03",
    "geos": "#B0B0B0",
    "ifpen": "#000000",
    "kfupm": "#FF7573",
    "opengosim": "#AD92C3",
    "opm": "#7141A3",
    "pau-inria": "#AC7A71",
    "pflotran": "#683229",
    "rice": "#F39ABF",
    "sintef": "#D651A9",
    "slb": "#75D166",
    "stuttgart": "#127F12",
    "tetratech": "#07A4BB",
    "ut-csee": "#A1A30C"
}

mass_cmap = matplotlib.colormaps['icefire'].resampled(256)
mass_cmap.__dict__['colors'][0] = [1, 1, 1, 1]
mass_cmap.set_bad([0.5, 0.5, 0.5])
mass_cmap.set_under([1, 1, 1])

def plotColors():
    fig, axs = plt.subplots(figsize=(6, 3))

    dx = 0.6
    pad = 0.5*(1 - dx)

    for i, group in zip(range(len(groups_and_colors)), groups_and_colors):
        color = groups_and_colors[group]
        x = (i%6) + pad
        y = int(i/6) + pad

        rect = matplotlib.patches.Rectangle((x, y), dx, dx, color=color)
        axs.add_patch(rect)
        axs.text(x + 0.5*dx, y - 1.1*pad, group, horizontalalignment='center', verticalalignment='bottom')

    axs.set_xlim((0, 6))
    axs.set_ylim((0, 3))
    axs.axis('off')

    fig.savefig('groups_and_colors.png', bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    plotColors()
