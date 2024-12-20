# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import seaborn

groups_and_colors = {
    "calgary": "#000000",
    "cau-kiel": "#A30059",
    "csiro": "#1B4400",
    "ctc-cne": "#1CE6FF",
    "darts": "#FF34FF",
    "geos": "#FF4A46",
    "ifpen": "#008941",
    "kfupm": "#006FA6",
    "opengosim": "#EFCBD5",#"#FFDBE5",
    "opm": "#7A4900",
    "pau-inria": "#0000A6",
    "pflotran": "#63FFAC",
    "rice": "#004D43",
    "sintef": "#8FB0FF",
    "slb": "#997D87",
    "stuttgart": "#5A0007",
    "tetratech": "#809693",
    "ut-csee": "#4FC601"
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
