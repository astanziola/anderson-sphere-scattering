%% Point Source Scattering Example
% This script demonstrates the acoustic scattering from a fluid sphere using
% Anderson's analytical solution. It visualizes both the real part and magnitude
% of the scattered pressure field in the xz-plane.

%% Clean workspace
clear
close all

%% Add functions to path
addpath(genpath("."))

%% Define physical parameters
% Medium properties (normalized units)
c0 = 1.0;        % Speed of sound in surrounding medium
rho0 = 1.0;      % Density of surrounding medium

% Sphere properties
c1 = 1.5;        % Speed of sound in sphere
rho1 = 3.0;      % Density of sphere
R = 0.35;        % Sphere radius

% Wave properties
f = 15;         % Base frequency
omega = 2*pi*f;  % Angular frequency
order = 250;     % Maximum order for modal expansion
D = 1.2;         % Point source distance on the z-axis

% Wavelength in the background medium
fprintf('Wavelength in the background medium: %.3f m\n', c0/f)

%% Setup computational grid
% Define spatial grid in xz-plane (y = 0)
grid_points = 301;
domain_size = 1.0;
x_range = linspace(-domain_size, domain_size, grid_points);
z_range = linspace(-domain_size, domain_size, grid_points);
[X, Z] = meshgrid(x_range, z_range);

%% Create position matrix for vectorized computation
% Reshape grid into a matrix of 3D points (y = 0)
positions = zeros(3, numel(X));
positions(1,:) = X(:)';  % x coordinates
positions(2,:) = 0;      % y coordinates (all zero)
positions(3,:) = Z(:)';  % z coordinates

%% Calculate pressure field (vectorized)
% Compute pressure for all valid points at once
P = computeAndersonSphereSolution(positions, c0, rho0, c1, rho1, R, omega, order, ...
    'kind', 'point-source', 'D', D);
P = reshape(P, size(X));

%% Setup visualization
figure('Position', [100 100 1500 400], 'Name', 'Anderson Sphere Scattering (Vectorized)')

% Setup real part subplot
ax(1) = subplot(1,3,1);
h1 = imagesc(x_range, z_range, real(P));
axis equal
colormap(ax(1),'jet')
cb1 = colorbar;
ylabel(cb1, 'Pressure (Real Part)')
title('Real Part of Pressure Field')
xlabel('x (m)')
ylabel('z (m)')
hold on

% Draw sphere outline
theta = linspace(0, 2*pi, 100);
plot(R*cos(theta), R*sin(theta), 'w.', 'LineWidth', 3)
hold off

% Setup magnitude subplot
ax(2) = subplot(1,3,2);
h2 = imagesc(x_range, z_range, abs(P));
axis equal
colormap(ax(2),'hot')
cb2 = colorbar;
ylabel(cb2, 'Pressure (Magnitude)');
title('Absolute Value of Pressure Field')
xlabel('x (m)')
ylabel('z (m)')
hold on

% Draw sphere outline
theta = linspace(0, 2*pi, 100);
plot(R*cos(theta), R*sin(theta), 'w.', 'LineWidth', 3)
hold off

% Setup z-axis line plot
ax(3) = subplot(1,3,3);
center_i = ceil(grid_points/2);
p_z_line = P(:, center_i);
plot(z_range, abs(p_z_line), 'b-', 'LineWidth', 1.5);
hold on
plot(z_range, real(p_z_line), 'r--', 'LineWidth', 1.5);
xlabel('z (m)')
ylabel('Pressure')
title('Pressure Along z-axis (x = 0)')
legend('|P|', 'Re(P)', 'Sphere', 'Location', 'best')
grid on

% Set color limits for consistent visualization
p_real_max = max(abs(real(P(:))));
p_abs_max = max(abs(P(:)));
p_abs_min = min(abs(P(:)));

subplot(1,3,1); clim([-p_real_max, p_real_max]);
subplot(1,3,2); clim([p_abs_min, p_abs_max]);
subplot(1,3,3); ylim([-p_abs_max, p_abs_max]);

%% Add title with simulation parameters
sgtitle(sprintf(['Acoustic Scattering from Fluid Sphere\n', ...
                 'c₁/c₀ = %.1f, ρ₁/ρ₀ = %.1f, kR = %.1f'], ...
                 c1/c0, rho1/rho0, omega/c0*R), ...
        'FontSize', 12);
