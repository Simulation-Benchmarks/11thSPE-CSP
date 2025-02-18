# SPDX-FileCopyrightText: 2024 Bernd Flemisch <bernd.flemisch@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
#!/usr/bin/env python3

def add_legend(ax):
    ax.plot([], [], color='k', linestyle='-', label='1')
    ax.plot([], [], color='k', linestyle='--', label='2')
    ax.plot([], [], color='k', linestyle='-.', label='3')
    ax.plot([], [], color='k', linestyle=':', label='4')

    numAxes = len(ax.get_figure().axes)
    handles, labels = ax.get_legend_handles_labels()

    by_label = dict(zip(labels[-4:], handles[-4:]))
    leg = ax.legend(by_label.values(), by_label.keys(), loc='upper left', title='result')
    ax.add_artist(leg)

    by_label = dict(zip(labels[:-4], handles[:-4]))
    if numAxes == 1:
        leg = ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1, 0.5))
    elif numAxes == 2:
        leg = ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(1.09, -0.2), ncols=8)
    elif numAxes == 4:
        leg = ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(1.09, -1.42), ncols=7)
    else:
        leg = ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(1.09, -2.6), ncols=10)
    for line in leg.get_lines():
        line.set_linestyle('-')
        line.set_linewidth(4)
        if line.get_label() == 'from dense data':
            line.set_linewidth(0)
            line.set_color('k')
            line.set_marker('$\u2736$')
