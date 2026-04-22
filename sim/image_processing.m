function [img_matrix, sumx, sumy] = image_processing(Cell_velocity, I2_list, ...
        rep_std, rep_vel, vel_pix, atten_factorx, atten_factory)
% IMAGE_PROCESSING  Apply recoil noise, PSF convolution, and generate projections.
%
% Electron recoil velocity is added before histogramming.
% The 2-D Voigt PSF is convolved with each histogram independently.
% parfor is used across intensity groups (each group is fully independent).
%
% Inputs:
%   Cell_velocity    — {n_groups×1} cell of MAX_per_group×3 velocity arrays [m/s]
%   I2_list          — n_groups×1 vector of I2 intensity values (for indexing)
%   rep_std          — recoil velocity std [m/s]
%   rep_vel          — recoil velocity magnitude [m/s]
%   vel_pix          — velocity per pixel [m/s/pixel]
%   atten_factorx    — PSF x-axis attenuation factor
%   atten_factory    — PSF y-axis attenuation factor
%
% Outputs:
%   img_matrix — {n_groups×1} cell of 251×251 blurred histogram images
%   sumx       — 251 × n_groups x-projection (sum over y)
%   sumy       — 251 × n_groups y-projection (sum over x)

    n_groups     = numel(I2_list);
    MAX_per_group = size(Cell_velocity{1}, 1);

    %% Electron recoil noise — same random kicks for all groups (correlated)
    vx_rand = rep_std * randn(MAX_per_group, 1);
    vy_rand = rep_std * randn(MAX_per_group, 1);
    vz_rand = rep_std * randn(MAX_per_group, 1);
    v_norm  = sqrt(vx_rand.^2 + vy_rand.^2 + vz_rand.^2);

    Cell_velocity_el = cell(n_groups, 1);
    for g = 1:n_groups
        Cell_velocity_el{g}(:,1) = Cell_velocity{g}(:,1) + rep_vel * vx_rand ./ v_norm;
        Cell_velocity_el{g}(:,2) = Cell_velocity{g}(:,2) + rep_vel * vy_rand ./ v_norm;
    end

    %% PSF (2-D Voigt) — computed once, reused for all groups
    xgrid = -125:125;   % pixel grid
    ygrid = -125:125;
    [xgrid_2, ygrid_2] = meshgrid(xgrid, ygrid);

    pFit = [-3.698853876412670,  2.554394107105213e+04, 0.441357801223283, ...
             6.354598584616713,  6.584074190809296,     8.250163296670557, ...
             9.417167979585598,  0,                     0, -0.477546278665701];
    M = voigt2d_aniso_rot(xgrid_2 * atten_factorx, ygrid_2 * atten_factory, pFit);

    %% Convert velocities to pixel units (done outside parfor for clarity)
    vx_pix = cell(n_groups, 1);
    vy_pix = cell(n_groups, 1);
    for g = 1:n_groups
        vx_pix{g} = Cell_velocity_el{g}(:,1) ./ vel_pix;
        vy_pix{g} = Cell_velocity_el{g}(:,2) ./ vel_pix;
    end

    %% Histogram + PSF convolution — each group is independent → parfor
    img_matrix = cell(n_groups, 1);

    parfor g = 1:n_groups
        h = hist3([vx_pix{g}, vy_pix{g}], 'Ctrs', {xgrid, ygrid});
        img_matrix{g} = conv2(h, M, 'same');
    end

    %% Projections (x and y marginals)
    sumx = zeros(numel(xgrid), n_groups);
    sumy = zeros(numel(ygrid), n_groups);
    for g = 1:n_groups
        sumx(:,g) = sum(img_matrix{g}, 2);
        sumy(:,g) = sum(img_matrix{g}, 1);
    end

end
