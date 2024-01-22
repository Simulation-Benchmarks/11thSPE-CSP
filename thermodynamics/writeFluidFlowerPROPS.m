% SPDX-FileCopyrightText: 2023 Olav Møyner <Olav.Moyner@sintef.no>
%
% SPDX-License-Identifier: MIT
function writeFluidFlowerPROPS(tab_h2o, tab_co2, tab_sol, varargin)
    % Write FluidFlower properties as black oil tables.
    % Written by Olav Møyner. Copyright SINTEF AS (2023).
    % SPDX-License-Identifier: MIT
    opts = struct(...
        'rs', true,        ... % CO2 can dissolve into liquid phase
        'rv', true,        ... % Water can vaporize into vapor phase
        'plot', false,     ... % Plot the tables for debugging
        'dir', '',         ... % Output folder to write files to
        'units', 'metric', ... % Unit system (si, lab, metric or field)
        'cap', true,       ... % Ensure that Rs -> 0 as p -> 0
        'nusat', 5         ... % Number of undersaturated points
    );
    for i = 1:(numel(varargin)/2)
        key = varargin{2*(i-1)+1};
        val = varargin{2*(i-1)+2};
        assert(isfield(opts, key));
        opts.(key) = val;
    end
    barsa = 1e5;
    u = unitFactors(opts.units);
    if ~isempty(opts.dir)
        if ~isfolder(opts.dir)
            mkdir(opts.dir);
        end
    end
    do_plot = opts.plot;
    nusat = opts.nusat;

    dissolve_gas = opts.rs;
    vaporize_water = opts.rv;
    T = tab_co2.x_Temperature__C_(1);
    % 
    tab_sol = tab_sol(tab_sol.x_Temperature__C_ == T, :);
    tab_h2o = tab_h2o(tab_h2o.x_Temperature__C_ == T, :);
    tab_co2 = tab_co2(tab_co2.x_Temperature__C_ == T, :);

    p_sol = tab_sol.phasePressure_Pa_;
    p_h2o = tab_h2o.pressure_Pa_;
    p_co2 = tab_co2.pressure_Pa_;

    assert(all(p_sol == p_h2o));
    assert(all(p_sol == p_co2));

    p = p_sol/barsa;
   
    molar_mass_co2 = 44.1e-3; % kg^3/mol
    molar_mass_h2o = 18.01528e-3; % kg^3/mol

    % Gas in liquid phase
    x_co2 = tab_sol.x_CO2___;
    x_co2 = x_co2.*dissolve_gas;
    x_co2 = max(x_co2, 1e-8);
    % Water in gas phase
    y_h2o = tab_sol.y_H2O___;
    y_h2o = y_h2o.*vaporize_water;

    X_co2 = mole_fraction_to_mass_fraction(x_co2, molar_mass_co2, molar_mass_h2o);
    Y_h2o = mole_fraction_to_mass_fraction(y_h2o, molar_mass_h2o, molar_mass_co2);

    if do_plot
        figure(1); clf;
        subplot(1, 2, 1)
        plot(p_sol, tab_sol.x_CO2___);
        title('x CO2');
        subplot(1, 2, 2)
        plot(p_sol, tab_sol.y_H2O___);
        title('y H2o');
    end
    % First plot the properties for the benchmark
    pure_water_density = tab_h2o.density_kg_m3_;
    pure_gas_density = tab_co2.density_kg_m3_;

    water_viscosity = tab_h2o.viscosity_Pa_s_;
    gas_viscosity = tab_co2.viscosity_Pa_s_;

    saturated_water_density = water_density(T, pure_water_density, X_co2);

    if do_plot

        figure(1); clf; hold on
        subplot(1, 2, 1);  hold on
        plot(p, pure_water_density);
        plot(p, saturated_water_density);
        legend('Pure', 'Saturated with CO2');

        title('Water density');
        subplot(1, 2, 2); hold on
        plot(p, pure_gas_density);
        title('Gas density');
    end
    %
    rhoGS = tab_co2.density_kg_m3_(1);
    rhoOS = tab_h2o.density_kg_m3_(1);

    % Make sure we don't use the tables beyond this point
    clear tab_h2o tab_co2 tab_sol

    % R_s sat = rhoOS * mass_fraction_co2_in_liquid_phase / 
    %          (rhoGS - mass_fraction_co2_in_liquid_phase)
    R_s = rhoOS.*X_co2./(rhoGS.*(1 - X_co2));
    % R_v sat = rhoGS * mass_fraction_h2o_in_vapor_phase /
    %          (rhoOS - mass_fraction_co2_in_vapor_phase)
    R_v = rhoGS.*Y_h2o./(rhoOS.*(1 - Y_h2o));
    if do_plot
        figure(1); clf; hold on
        subplot(1, 2, 1)
        plot(p_sol, R_s);
        title('Saturated R_s');
        subplot(1, 2, 2)
        plot(p_sol, R_v);
        title('Saturated R_v');
    end

    % Convert <-> mass fraction to Rs/Rv using that formula
    % From that, calculate densities etc

    % b-factors: bO(p, Rs) = rhoO(p, T, X)/(rhoOS + Rs(X)*rhoGS)
    % b-factors: bG(p, Rv) = rhoG(p, T, X)/(rhoGS + Rv(X)*rhoOS)

    bG = shrinkage_factor(pure_gas_density, Y_h2o, rhoGS, rhoOS);
    bO = shrinkage_factor(water_density(T, pure_water_density, X_co2), Y_h2o, rhoOS, rhoGS);

    if do_plot
        figure(1); clf; hold on
        subplot(1, 2, 1)
        plot(p_sol, bG);
        title('Saturated b_g');
        subplot(1, 2, 2)
        plot(p_sol, bO);
        title('Saturated b_o');
    end

    if dissolve_gas
        write_pvto(p_sol, R_s, pure_water_density, water_viscosity, T, X_co2, Y_h2o, rhoOS, rhoGS, u, nusat, opts)
    else
        write_immiscible(p_sol, pure_water_density./rhoOS, water_viscosity, u, 'PVDO', opts);
    end
    if vaporize_water
        write_pvtg(p_sol, R_v, pure_gas_density, gas_viscosity, Y_h2o, rhoGS, rhoOS, u, opts)
    else
        write_immiscible(p_sol, pure_gas_density./rhoGS, gas_viscosity, u, 'PVDG', opts);
    end
    file_path = fullfile(opts.dir, 'DENSITY.txt');
    fn = fopen(file_path, 'w');
    fprintf(fn, 'DENSITY\n    %f 1 %f /\n', rhoOS, rhoGS);
    fprintf('%s written.\n\n', file_path);
end

function write_pvto(p, Rs, pure_water_density, water_viscosity, T, X_co2_sat, Y_h2o_sat, rhoOS, rhoGS, u, n, opts)
    file_path = fullfile(opts.dir, 'PVTO.txt');
    fn = fopen(file_path, 'w');
    % PVTO: Tables by R_s: For each, given a set of p_o B_o mu_o
    fprintf(fn, '-- RS    PRES    FVF      VISC\n');
    fprintf(fn, 'PVTO\n');

    if opts.cap
        X_co2_min = 1e-8;
        Rs_min = rhoOS.*X_co2_min./rhoGS;
        p_zero = interp1(X_co2_sat, p, X_co2_min, 'linear', 'extrap');
        assert(p_zero > 0)

        rhoW_min = interp1(p, pure_water_density, p_zero, 'linear', 'extrap');
        p = [p_zero; p];
        Rs = [Rs_min; Rs];
        X_co2_sat = [X_co2_min; X_co2_sat];
        Y_h2o_sat = [Y_h2o_sat(1); Y_h2o_sat]; % Note: May be wrong for Rs+Rv case...
        pure_water_density = [rhoW_min; pure_water_density];
        water_viscosity = [water_viscosity(1); water_viscosity];
    end
    bO_sat = shrinkage_factor(water_density(T, pure_water_density, X_co2_sat), X_co2_sat, rhoOS, rhoGS);

    M = numel(Rs);
    uk = u.gasvol_s / u.liqvol_s;
    ub = u.liqvol_r/u.liqvol_s;
    for i = 1:M
        mu = water_viscosity(i);
        fprintf(fn, '%-4.10g %-4.10e %-4.10g %-4.10g\n', Rs(i)*uk, p(i)*u.press, ub*(1/bO_sat(i)), u.viscosity*mu);
        X_co2 = X_co2_sat(i);
        % Y_h2o = Y_h2o_sat(i);
        p_usat = linspace(p(i), 1.2*max(p), n+1);
        p_usat = p_usat(2:end);
        for j = 1:n
            rho_w_j = interp1(p, pure_water_density, p_usat(j), 'linear', 'extrap');
            rho_w_dissolved = water_density(T, rho_w_j, X_co2);
            bO_usat = shrinkage_factor(rho_w_dissolved, X_co2, rhoOS, rhoGS);
            rho_new = bO_usat*(rhoOS + Rs(i)*rhoGS);
            assert(abs(rho_w_dissolved - rho_new) < 1e-3)
            fprintf(fn, '           %-4.10g %-4.10g %-4.10g', p_usat(j)*u.press, ub*(1/bO_usat), u.viscosity*mu);
            if j == n
                fprintf(fn, ' /\n');
            else
                fprintf(fn, '\n');
            end
        end
    end
    fprintf(fn, '/\n');
    fclose(fn);
    fprintf('%s written.\n', file_path);
end

function write_pvtg(p, Rv, pure_gas_density, gas_viscosity, Y_h2o, rhoGS, rhoOS, u, opts)
    bG = shrinkage_factor(pure_gas_density, Y_h2o, rhoGS, rhoOS);
    file_path = fullfile(opts.dir, 'PVTG.txt');
    fn = fopen(file_path, 'w');

    uk = u.press;
    urv = u.liqvol_r/u.gasvol_s;
    ub = u.gasvol_r/u.gasvol_s;
    umu = u.viscosity;

    % PVTG: Tables by pressure: For each, given a set of R_v B_g mu_g
    fprintf(fn, '-- PRES  RV    FVF      VISC\n');
    fprintf(fn, 'PVTG\n');

    M = numel(Rv);
    for i = 1:M
        rv = Rv(i);
        mu = gas_viscosity(i);
        fprintf(fn, '%-4.10g %-4.10g %-4.10g %-4.10g\n', uk*p(i), urv*rv, ub/bG(i), umu*mu);
        fprintf(fn, '             %-4.10g %-4.10g %-4.10g', 0.0, ub/bG(i), umu*mu);
        fprintf(fn, ' /\n');
    end
    fprintf(fn, '/\n');
    fprintf('%s written.\n', file_path);
    fclose(fn);
end

function write_immiscible(p, b, mu, u, title, opts)
    filename = [title, '.txt'];
    file_path = fullfile(opts.dir, filename);
    fn = fopen(file_path, 'w');
    if strcmpi(title, 'pvdg')
        fvf_u = u.gasvol_r/u.gasvol_s;
    else
        fvf_u = u.liqvol_r/u.liqvol_s;
    end

    fprintf(fn, '-- PRES    FVF      VISC\n');
    fprintf(fn, '%s\n', title);
    for i = 1:numel(p)
        fprintf(fn, '%-4.10g %-4.10g %-4.10g\n', p(i)*u.press, fvf_u/b(i), u.viscosity*mu(i));
    end
    fprintf(fn, '/\n');
    fclose(fn);
    fprintf('%s written.\n', file_path);
end

function bO = shrinkage_factor(rhoO, X_co2, rhoOS, rhoGS)
    Rs = rhoOS.*X_co2./((1-X_co2).*rhoGS);
    bO = rhoO./(rhoOS + Rs*rhoGS);
end

function rho = co2_density_in_water(T)
    V = co2_partial_volume_in_water(T);
    M = 44.1e-3; % kg / mol
    rho = M/V;
end

function V = co2_partial_volume_in_water(T)
    % Simple check to avoid Kelvins
    assert(all(T < 273.15))
    V = (37.51 - T*9.585e-2 + (T^2)*8.74e-4 + (T^3)*5.044e-7)/1e6;
end

function rho = water_density(T, pure_water_density, X_co2)
    % M_co2_in_water / (M_co2_in_water + M_h2o_in_water) 
    X_h2o = 1 - X_co2;
    vol = X_h2o./pure_water_density + X_co2./co2_density_in_water(T);
    % rho_inv
    % pure_water_density
    %co2_density_in_water(T)
    rho = 1./vol;
end

function X = mole_fraction_to_mass_fraction(x, mw_self, mw_other)
    mass_self = x.*mw_self;
    mass_other = (1-x).*mw_other;
    
    X = mass_self./(mass_self + mass_other);
end

function u = unitFactors(units)
    % Simple version of unitConversionFactors from MRST.
    [mu, liqvol_s, gasvol_s, liqvol_r, gasvol_r, p] = deal(1);
    switch lower(units)
        case 'si'
            % Nothing to do.
        case 'field'
            mu = 1e-3; % cP
            p = 6.894757e3; % psi
            liqvol_s = 1.589872949280001e-01; % stb
            liqvol_r = liqvol_s;
            gasvol_s = 2.831684659200000e+01;
            gasvol_r = liqvol_s;
        case 'lab'
            p = 101325;
            liqvol_s = 1e-6;
            gasvol_s = 1e-6;
            liqvol_r = 1e-6;
            gasvol_r = 1-6;
            mu = 1e-3; % cP
        case 'metric'
            mu = 1e-3; % cP
            p = 1e5; % bar
        otherwise
            error('Unknown unit system')
    end
    u = struct('press', 1/p,      ...
               'viscosity', 1/mu, ...
               'liqvol_r', 1/liqvol_r,  ...
               'gasvol_r', 1/gasvol_r,  ...
               'liqvol_s', 1/liqvol_s,  ...
               'gasvol_s', 1/gasvol_s   ...
       );
end
