%% Anderson Sphere Scattering Convergence Study
% This script demonstrates how the scattered field converges along different
% lines (z-axis, x-axis, and diagonal) as the order of modal expansion increases.

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
c1 = 1.7;        % Speed of sound in sphere
rho1 = 2.5;      % Density of sphere
R = 0.2;         % Sphere radius

% Wave properties
f = 30.0;         % Base frequency
omega = f*2*pi; % Angular frequency (chosen for visible scattering effects)
max_order = 100; % Maximum order for convergence study

%% Setup computational grid
% Define spatial points
grid_points = 401;
domain_size = 0.5;
line_range = linspace(-domain_size, domain_size, grid_points);
orders = 0:max_order;

% Initialize pressure field matrices for each line
P_z = zeros(length(line_range), length(orders));  % z-axis (x=y=0)
P_x = zeros(length(line_range), length(orders));  % x-axis (y=z=0)
P_diag = zeros(length(line_range), length(orders));  % diagonal (x=z, y=0)

%% Setup visualization
figure('Position', [100 100 1500 500], 'Name', 'Anderson Sphere Convergence Study')

% Create subplots for each line
subplot_titles = {'Z-axis (x=y=0)', 'X-axis (y=z=0)', 'Diagonal (x=z, y=0)'};
ax = zeros(1,3);
h = zeros(1,3);
cb = zeros(1,3);

for i = 1:3
    ax(i) = subplot(1,3,i);
    h(i) = imagesc(orders, line_range, zeros(length(line_range), length(orders)));
    colormap(ax(i), 'turbo')
    cb(i) = colorbar;
    title(cb(i), 'Pressure Magnitude')
    xlabel('Order of Modal Expansion')
    ylabel('Position (m)')
    title(subplot_titles{i})
    
    % Add sphere intersection markers
    hold(ax(i), 'on')
    if i == 1 || i == 2  % For z-axis and x-axis
        plot(ax(i), [orders(1) orders(end)], [R R], 'w--', 'LineWidth', 1.5);
        plot(ax(i), [orders(1) orders(end)], [-R -R], 'w--', 'LineWidth', 1.5);
    else  % For diagonal
        R_diag = R/sqrt(2);  % Intersection points along diagonal
        plot(ax(i), [orders(1) orders(end)], [R_diag R_diag], 'w--', 'LineWidth', 1.5);
        plot(ax(i), [orders(1) orders(end)], [-R_diag -R_diag], 'w--', 'LineWidth', 1.5);
    end
    hold(ax(i), 'off')
end

%% Calculate pressure field evolution
for n = 1:length(orders)
    current_order = orders(n);
    
    % Calculate pressure along each line for current order
    for i = 1:length(line_range)
        pos = line_range(i);
        
        % Z-axis points (x=y=0)
        if abs(pos) > R
            p = [0; 0; pos];
            P_z(i,n) = computeAndersonSphereSolution(p, c0, rho0, c1, rho1, R, omega, current_order);
        end
        
        % X-axis points (y=z=0)
        if abs(pos) > R
            p = [pos; 0; 0];
            P_x(i,n) = computeAndersonSphereSolution(p, c0, rho0, c1, rho1, R, omega, current_order);
        end
        
        % Diagonal points (x=z, y=0)
        if abs(pos) > R/sqrt(2)  % Account for diagonal intersection
            p = [pos; 0; pos];  % x = z
            P_diag(i,n) = computeAndersonSphereSolution(p, c0, rho0, c1, rho1, R, omega, current_order);
        end
    end
    
    % Update all visualizations
    set(h(1), 'CData', abs(P_z));
    set(h(2), 'CData', abs(P_x));
    set(h(3), 'CData', abs(P_diag));
    

    drawnow limitrate
end

% Final update
txt.String = 'Calculation Complete';

%% Add title with simulation parameters
sgtitle(sprintf(['Convergence of Pressure Field Along Different Lines\n', ...
                 'c₁/c₀ = %.1f, ρ₁/ρ₀ = %.1f, kR = %.1f'], ...
                 c1/c0, rho1/rho0, omega/c0*R), ...
         'FontSize', 12)

%% Add convergence analysis plot
figure('Position', [100 600 1500 400], 'Name', 'Convergence Analysis')

% Calculate relative change between successive orders for each line
rel_change_z = zeros(1, length(orders)-1);
rel_change_x = zeros(1, length(orders)-1);
rel_change_diag = zeros(1, length(orders)-1);

for n = 1:length(orders)-1
    rel_change_z(n) = norm(abs(P_z(:,n+1) - P_z(:,n))) / norm(abs(P_z(:,n)));
    rel_change_x(n) = norm(abs(P_x(:,n+1) - P_x(:,n))) / norm(abs(P_x(:,n)));
    rel_change_diag(n) = norm(abs(P_diag(:,n+1) - P_diag(:,n))) / norm(abs(P_diag(:,n)));
end

% Plot convergence
semilogy(orders(2:end), rel_change_z, 'b-', 'LineWidth', 1.5)
hold on
semilogy(orders(2:end), rel_change_x, 'r-', 'LineWidth', 1.5)
semilogy(orders(2:end), rel_change_diag, 'g-', 'LineWidth', 1.5)
yline(1e-6, 'k--', 'Typical convergence threshold', 'LineWidth', 1.5)
hold off

grid on
xlabel('Order of Modal Expansion')
ylabel('Relative Change')
title('Convergence Analysis')
legend('Z-axis', 'X-axis', 'Diagonal', 'Location', 'southwest')