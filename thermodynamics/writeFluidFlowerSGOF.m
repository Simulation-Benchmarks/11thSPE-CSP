function writeFluidFlowerSGOF(varargin)
    % Write fluid flower SGOF table.
    % Written by Olav Møyner. Copyright SINTEF AS (2023).
    % SPDX-License-Identifier: MIT
    opt = struct('n', 100, ...
                 'plot', false, ...
                 'dir', '', ...
                 'filename', 'SGOF.txt', ...
                 'units', 'metric');
    for i = 1:(numel(varargin)/2)
        key = varargin{2*(i-1)+1};
        val = varargin{2*(i-1)+2};
        assert(isfield(opts, key));
        opts.(key) = val;
    end
    u = unitFactors(opt.units);

    if ~isempty(opt.dir)
        if ~isfolder(opt.dir)
            mkdir(opt.dir);
        end
    end
    nreg = 6;
    sw_imm = [0.32, 0.14, 0.12, 0.12, 0.12, 0.10];
    sg_imm = repmat(0.1, 1, 6);
    c_a1 = 2;
    c_a2 = 2;
    Pe = [1500, 300, 100, 25, 10, 1];
    P_c_max = repmat(9.5e4, 1, 6);
    
    if opt.plot
        figure(1);
        clf;
        figure(2);
        clf;
    end
    fn = fopen(fullfile(opt.dir, opt.filename), 'w');
    fprintf(fn, '-- SG    KRG    KROG   PCOG\n');
    fprintf(fn, 'SGOF\n');

    for reg = 1:nreg
        s0 = sg_imm(reg);
        sat = [0, linspace(s0, 1, opt.n-1)];
        sg = scaleSaturation(sat, s0);
        so = scaleSaturation(1-sat, sw_imm(reg));
        krg = sg.^c_a1;
        kro = so.^c_a1;
        pc_og_bar = Pe(reg)*so.^(1-c_a2);
        
        pcmax = P_c_max(reg);
        pc_og = pcmax*erf((pc_og_bar/pcmax)*(sqrt(pi)/2));
        if opt.plot
            figure(1);
            subplot(1, 2, 1); hold on;
            title('krO')
            plot(1 - sat, kro);
            subplot(1, 2, 2); hold on;
            title('krG')
            plot(sat, krg);
            figure(2);
            subplot(1, nreg, reg); hold on;
            title('pcOG')
            plot(sat, pc_og_bar, 'o');
            plot(sat, pc_og);
        end
        for i = 1:opt.n
            fprintf(fn, '%1.8f %1.8f %1.8f %1.8f\n', sat(i), krg(i), kro(i), pc_og(i)*u.press);
        end
        fprintf(fn, '/\n');
    end
    disp([opt.filename ' written.']);
    fclose(fn);
end

function s_scaled = scaleSaturation(s, s_imm)
    s_scaled = max((s - s_imm)/(1 - s_imm), 0);
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
