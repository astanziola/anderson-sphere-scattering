%% Anderson Sphere Scattering Example
% This script demonstrates the acoustic scattering from a fluid sphere using
% Anderson's analytical solution. It visualizes both the real part and magnitude
% of the scattered pressure field in the xz-plane.
%
% The visualization is created progressively, starting from the center and working
% outward, to provide a real-time view of the computation progress.

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
c1 = 2.0;        % Speed of sound in sphere
rho1 = 2.0;      % Density of sphere
R = 0.25;         % Sphere radius

% Wave properties
f = 5.0;         % Base frequency
omega = 2*pi*f;   % Angular frequency
order = 100;       % Maximum order for modal expansion
D = 1.2;          % Point source distance on the z-axis

% Wavelength in the background medium
fprintf('Wavelength in the background medium: %.3f m\n', c0/f)

%% Setup computational grid
% Define spatial grid in xz-plane (y = 0)
grid_points = 201;
domain_size = 1.0;
x_range = linspace(-domain_size, domain_size, grid_points);
z_range = linspace(-domain_size, domain_size, grid_points);
[X, Z] = meshgrid(x_range, z_range);

% Initialize pressure field
P = zeros(size(X));

%% Setup visualization
figure('Position', [100 100 1500 400], 'Name', 'Anderson Sphere Scattering')

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
plot(R*cos(theta), R*sin(theta), 'w.', 'LineWidth', 3)
hold off

% Setup z-axis line plot
ax(3) = subplot(1,3,3);
h3 = plot(z_range, zeros(size(z_range)), 'b-', 'LineWidth', 1.5);
hold on
plot(z_range, zeros(size(z_range)), 'r--', 'LineWidth', 1.5);
plot([-R R], [0 0], 'k', 'LineWidth', 20);  % Sphere boundaries
xlabel('z (m)')
ylabel('Pressure')
title('Pressure Along z-axis (x = 0)')
legend('|P|', 'Re(P)', 'Sphere', 'Location', 'best')
grid on
hold off

%% Calculate pressure field
% Process points in order of increasing distance from center for better visualization
center_i = ceil(grid_points/2);
center_j = ceil(grid_points/2);

% Create flat indices and calculate distances from center
[X_flat, Z_flat] = meshgrid(1:grid_points, 1:grid_points);
X_flat = X_flat(:);
Z_flat = Z_flat(:);
distances = sqrt((X_flat - center_i).^2 + (Z_flat - center_j).^2);
[~, order_indices] = sort(distances);

% Setup progress display
num_points = numel(X);
point_count = 0;
txt = annotation('textbox', [0.4, 0.95, 0.2, 0.05], ...
                'String', 'Progress: 0%', ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center');

% Process points and update visualization
update_interval = 10;  % Update display every N points
p_abs_min = inf;
for idx = 1:length(order_indices)
    i = X_flat(order_indices(idx));
    j = Z_flat(order_indices(idx));
    
    % Calculate pressure only for points outside the sphere
    if sqrt(X(j,i)^2 + Z(j,i)^2) > R
        % Calculate pressure at current point
        p = [X(j,i); 0; Z(j,i)];  % Point in 3D space (y = 0)

        % Main computation
        P(j,i) = computeAndersonSphereSolution(p, c0, rho0, c1, rho1, R, omega, order, ...
            kind='point-source', D=D);

        % Update visualization periodically
        if mod(point_count, update_interval) == 0
            % Update pressure field plots
            set(h1, 'CData', real(P));
            set(h2, 'CData', abs(P));
            
            % Update color limits for consistent visualization
            p_real_max = max(abs(real(P(:))));
            p_abs_max = max(abs(P(:)));
            if abs(P(j,i)) < p_abs_min
                p_abs_min = abs(P(j,i));
            end

            if ~any(isnan([p_real_max, p_abs_max]))
                subplot(1,3,1); clim([-p_real_max, p_real_max]);
                if p_abs_max > p_abs_min
                    subplot(1,3,2); clim([p_abs_min, p_abs_max]);
                end
                
                % Update z-axis line plot
                z_line_idx = center_i * ones(size(z_range));
                p_z_line = P(:, center_i);
                set(h3, 'XData', z_range, 'YData', abs(p_z_line));
                h3.Parent.Children(2).YData = real(p_z_line);
                h3.Parent.YLim = [-p_real_max, p_real_max];
            end
            
            % Update progress display
            progress = point_count / num_points * 100;
            txt.String = sprintf('Progress: %.1f%%', progress);
            drawnow limitrate
        end
    end
    
    point_count = point_count + 1;
end

% Final update
txt.String = '';
drawnow

%% Add title with simulation parameters
sgtitle(sprintf(['Acoustic Scattering from Fluid Sphere\n', ...
                 'c₁/c₀ = %.1f, ρ₁/ρ₀ = %.1f, kR = %.1f'], ...
                 c1/c0, rho1/rho0, omega/c0*R), ...
        'FontSize', 12);