% SPDX-FileCopyrightText: 2023 Olav Møyner <Olav.Moyner@sintef.no>
%
% SPDX-License-Identifier: MIT
%
% This file should be run from the folder where you store your outputs from
% the python scripts (make_component_table / make_solubility_table).
Atm = 101325.0;
Bar = 1e5;
np = 50; % Number of pressure points
p_min = 1e-3*Atm; % Tiny value to make sure that the table is safe to extrapolate a bit
p_max = 2.5*Bar;
T = 20; % C

clc
fprintf('To generate input tables, run the following:\n');
for i = 1:2
    if i == 1
        name = '-cH2O';
    else
        name = '-cCO2';
    end
    fprintf('python ./make_component_table.py -t1 %f -t2 %f -nt 2 -p1 %f -p2 %f -np %d %s\n', T, 1.1*T, p_min, p_max, np, name)
end
fprintf('python ./make_solubility_table.py -t1 %f -t2 %f -nt 2 -p1 %f -p2 %f -np %d\n\n', T, 1.1*T, p_min, p_max, np)
warning('off')
tab_h2o = readtable('h2ovalues.csv');
tab_co2 = readtable('co2values.csv');
tab_sol = readtable('solubilities.csv');
warning('on')
disp('Writing mutally miscible table (both pseudocomponents in both phases)');
writeFluidFlowerPROPS(tab_h2o, tab_co2, tab_sol, 'rs', true, 'rv', true, 'dir', 'both_miscible');
disp('Writing one-way miscible table (CO2 dissolves into liquid, brine exists only as liquid)');
writeFluidFlowerPROPS(tab_h2o, tab_co2, tab_sol, 'rs', true, 'rv', false, 'dir', 'co2_miscible');
disp('Writing immiscible PVT');
writeFluidFlowerPROPS(tab_h2o, tab_co2, tab_sol, 'rs', false, 'rv', false, 'dir', 'immiscible');
disp('Writing saturation functions.');
writeFluidFlowerSGOF()
