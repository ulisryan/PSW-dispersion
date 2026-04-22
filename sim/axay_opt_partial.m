function [cos2theta, last_data, acceleration, temp] = axay_opt_partial( ...
        position, t, I1, I2, aligned_value, last_data, IR_position_stability, ...
        lambda, w0z, w0y, tau, MAX)
% AXAY_OPT_PARTIAL  Compute optical dipole force and acceleration for all particles.
%
% Inputs:
%   position              — MAX×3 [x, y, z] in meters
%   t                     — MAX×1 time per particle [s] (includes Ij jitter)
%   I1, I2                — MAX×1 (or scalar) peak laser intensities [W/m²]
%   aligned_value         — MAX×2001 cos2theta lookup table per particle
%   last_data             — MAX×2 [alpha_prev, I*gauss_prev] from previous step
%   IR_position_stability — MAX×2 [z_jitter, y_jitter] in meters
%   lambda                — laser wavelength [m]
%   w0z, w0y              — beam waist radii [m]
%   tau                   — laser pulse width [s]
%   MAX                   — number of particles
%
% Outputs:
%   cos2theta    — MAX×1 interpolated alignment values
%   last_data    — MAX×2 updated [alpha, I*gauss_t] for next step
%   acceleration — MAX×3 [ax, ay, az] in m/s²
%   temp         — MAX×7 debug array (matches original output format)

    % Physical constants
    epsilon0 = 8.854187e-12;        % vacuum permittivity [F/m]
    c        = 2.998e8;             % speed of light [m/s]
    inv_mass = 6.022e23 / 76.14e-3; % 1/m_CS2 [kg^-1]

    % Precompute temporal Gaussian envelope — evaluated once, reused 4× below
    gauss_t = exp(-4 * log(2) * t.^2 / tau^2);  % MAX×1

    % Spatial wave vector (x-direction fringe phase)
    q = 2 * pi * position(:,1) / lambda;

    % Interference fringe pattern (spatial, no time envelope)
    I_fringe  = I1 + I2 + 2 .* sqrt(I1 .* I2) .* cos(2.*q);

    % Gaussian beam profile (y and z transverse directions)
    spatial_z = exp(-2 .* (position(:,3) + IR_position_stability(:,1)).^2 / w0z^2);
    spatial_y = exp(-2 .* (position(:,2) + IR_position_stability(:,2)).^2 / w0y^2);

    % Total spatial intensity (no time envelope)
    I = I_fringe .* spatial_z .* spatial_y;

    % Scaled intensity for lookup table (units: 1e13 W/m²)
    I_scaled = I .* gauss_t .* 1e-13;
    I_floor  = floor(I_scaled);

    % Clamp column index to [1, n_cols-1] so interpolation (key+MAX) stays in bounds
    n_cols  = size(aligned_value, 2);
    col_idx = min(I_floor + 1, n_cols - 1);
    col_idx = max(col_idx, 1);

    % Linear interpolation of cos2theta from lookup table
    row_idx = (1:MAX)';
    key     = sub2ind(size(aligned_value), row_idx, col_idx);
    cos2theta = aligned_value(key) + ...
        (I_scaled - I_floor) .* (aligned_value(key + MAX) - aligned_value(key));

    % Molecular polarizability: linear function of cos2theta
    alpha = (16.8 - 6.2) * 1e-40 * cos2theta + 6.2e-40;

    % Radiation pressure prefactor: I / (2*epsilon0*c)
    mm = 0.5 * I / epsilon0 / c;

    % Dynamic polarizability response: d(alpha)/dI via finite difference from last step
    valid_prev = last_data(:,2) ~= 0;
    d_alpha_dI = (last_data(:,1) - alpha) ./ (last_data(:,2) - I .* gauss_t);
    d_alpha_dI = d_alpha_dI .* valid_prev;  % zero if no previous data available

    % Effective polarizability (static + dynamic response term)
    alpha_eff = alpha + I .* d_alpha_dI .* gauss_t;

    % Acceleration — only computed where I ~= 0 (GPU-compatible allocation via 'like')
    mask         = (I ~= 0);
    acceleration = zeros(MAX, 3, 'like', position);
    n_k          = 2 * pi / lambda;   % wavenumber

    % x: gradient force from interference fringe pattern
    x_grad = 2 .* n_k .* (2 .* sqrt(I1 .* I2) .* sin(2.*q)) ./ I_fringe;
    acceleration(mask, 1) = -mm(mask) .* x_grad(mask) .* inv_mass ...
                            .* gauss_t(mask) .* alpha_eff(mask);

    % y: gradient force from Gaussian beam waist (y-direction)
    acceleration(mask, 2) = -4 .* mm(mask) .* position(mask,2) .* inv_mass ...
                            .* gauss_t(mask) ./ w0y^2 .* alpha_eff(mask);

    % z: gradient force from Gaussian beam waist (z-direction)
    acceleration(mask, 3) = -4 .* mm(mask) .* position(mask,3) .* inv_mass ...
                            .* gauss_t(mask) ./ w0z^2 .* alpha_eff(mask);

    % Update state for next time step
    last_data(:,1) = alpha;
    last_data(:,2) = I .* gauss_t;

    % Debug output (format matches original)
    temp = [I, I_scaled, I_floor, row_idx, col_idx, key, cos2theta];
end
