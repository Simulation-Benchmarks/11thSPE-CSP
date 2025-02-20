# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3
import matplotlib.pyplot as plt

def set_fonts(size = 10):
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Computer Modern",
        "font.size": size,
        "legend.title_fontsize": "small",
        "legend.fontsize": "small"
    })

def add_legend(ax):
    for a in ax.get_figure().axes:
        if a.get_yaxis().get_scale() != 'log':
            a.ticklabel_format(style='sci', axis='y', scilimits=(-1,2))

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    numEntries = len(by_label)
    numAxes = len(ax.get_figure().axes)

    if numEntries%2 == 1 and numAxes > 1:
        ax.plot([], [], ' ', label=' ')
    if numEntries > 17 and numAxes == 4:
        ax.plot([], [], ' ', label='  ')
        ax.plot([], [], ' ', label='   ')
        ax.plot([], [], ' ', label='    ')
        ax.plot([], [], ' ', label='     ')
    if numEntries > 17 and numAxes == 2:
        ax.plot([], [], ' ', label='  ')
        ax.plot([], [], ' ', label='   ')
        ax.plot([], [], ' ', label='    ')
        ax.plot([], [], ' ', label='     ')
        ax.plot([], [], ' ', label='      ')
        ax.plot([], [], ' ', label='       ')
    ax.plot([], [], color='k', linestyle='-', label='result 1')
    ax.plot([], [], color='k', linestyle='--', label='result 2')
    ax.plot([], [], color='k', linestyle='-.', label='result 3')
    ax.plot([], [], color='k', linestyle=':', label='result 4')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    numEntries = len(by_label)

    if numAxes == 1:
        leg = ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1, 0.5))
    else:
        if numEntries < 7: ncols = 3
        elif numEntries < 9: ncols = 4
        elif numEntries < 11: ncols = 5
        elif numEntries < 13: ncols = 6
        else: ncols = 7
        if numAxes > 4 and numEntries > 14:
            if numEntries < 17: ncols = 8
            elif numEntries < 19: ncols = 9
            else: ncols = 10

        if numAxes == 2:
            leg = ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(1.11, -0.2), ncols=ncols)
        elif numAxes == 4:
            leg = ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(1.11, -1.42), ncols=ncols)
        else:
            leg = ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(1.11, -2.6), ncols=ncols)

    for line in leg.get_lines():
        if line.get_label() == 'from dense data':
            line.set_linewidth(0)
            line.set_color('k')
            line.set_marker('$*$')
        elif line.get_label()[:-2] != 'result' and line.get_label()[0] != ' ':
            line.set_linestyle('-')
            line.set_linewidth(4)
