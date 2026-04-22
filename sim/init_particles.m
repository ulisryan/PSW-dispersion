function [particle, laser] = init_particles(jm_data, MAX, ...
        x_position_fwhm, Dye_sigma, Dye_tau, ...
        vx_caused_tilted, x_init_vel_sigma, y_init_vel_sigma, ...
        z_init_vel_off,   z_init_vel_sigma, ...
        Ij_distribution,  Ips_z, Ips_y)
% INIT_PARTICLES  Initialize particle quantum states, positions, and velocities.
%
% Inputs:
%   jm_data           — path to temperature-specific alpha data file
%   MAX               — total number of particles (all groups combined)
%   x_position_fwhm  — initial x-position FWHM [m]
%   Dye_sigma         — Dye laser spatial jitter std [m]
%   Dye_tau           — Dye laser temporal jitter std [s]
%   vx_caused_tilted  — x-velocity offset from beam tilt [m/s]
%   x_init_vel_sigma  — x-velocity std [m/s]
%   y_init_vel_sigma  — y-velocity std [m/s]
%   z_init_vel_off    — z-velocity offset [m/s]
%   z_init_vel_sigma  — z-velocity std [m/s]
%   Ij_distribution   — IR timing jitter half-range [ps]
%   Ips_z, Ips_y      — IR position stability std [nm]
%
% Outputs:
%   particle — MAX×11 [x,y,z, vx,vy,vz, 0,0,0, j,|m|]
%   laser    — MAX×5  [y_jitter, t_jitter, Ij, Ips_z, Ips_y]

    %% Load quantum state distribution from temperature-specific file
    fid        = fopen(jm_data, 'r');
    alpha_data = fscanf(fid, '%lf', [4, Inf]);
    fclose(fid);
    alpha_data = alpha_data';

    input_j = alpha_data(:,1);
    input_m = alpha_data(:,2);
    proba   = alpha_data(:,4);

    % Convert cumulative probability to discrete probability
    proba(2:end) = proba(2:end) - proba(1:end-1);

    %% Sample quantum states weighted by Boltzmann distribution
    n_states     = numel(proba);
    random_index = randsample(1:n_states, MAX, true, proba);

    particle = zeros(MAX, 11);
    particle(:,10) = input_j(random_index);
    particle(:,11) = abs(input_m(random_index));

    %% Laser jitter and timing noise
    laser = zeros(MAX, 5);
    laser(:,1) = Dye_sigma * randn(MAX, 1);                                    % y spatial jitter [m]
    laser(:,2) = Dye_tau   * randn(MAX, 1) / sqrt(3);                          % temporal jitter [s]
    laser(:,3) = (-Ij_distribution + randi([0, Ij_distribution*2], MAX, 1)) * 1e-12;  % IR timing [s]
    laser(:,4) = Ips_z .* randn(MAX, 1) * 1e-9;                               % IR z-position stability [m]
    laser(:,5) = Ips_y .* randn(MAX, 1) * 1e-9;                               % IR y-position stability [m]

    %% Initial positions
    particle(:,1) = x_position_fwhm / 2.355 * randn(MAX, 1);  % x [m]
    particle(:,2) = Dye_sigma * randn(MAX, 1);                  % y [m]

    %% Initial velocities
    particle(:,4) = vx_caused_tilted + x_init_vel_sigma * randn(MAX, 1);  % vx [m/s]
    particle(:,5) = y_init_vel_sigma  * randn(MAX, 1);                     % vy [m/s]
    particle(:,6) = z_init_vel_off    + z_init_vel_sigma * randn(MAX, 1); % vz [m/s]

    %% Initial z-position (derived from velocity and timing jitter)
    particle(:,3) = (30e-9 + laser(:,2)) .* z_init_vel_off ...
                  - (60e-9 + laser(:,2)) .* particle(:,6) ...
                  + laser(:,1);

end
