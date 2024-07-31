clc
clear all
close all

Points = h5read("Structured_Mesh.h5", "/Points")';
Points_h = h5read("Structured_Mesh.h5", "/Points_h")';
Elements = h5read("Structured_Mesh.h5", "/Elements")' + 1;
Elements_h = h5read("Structured_Mesh.h5", "/Elements_h")' + 1;
JBP = h5read("Structured_Mesh.h5", "/JBP")';
JBV = h5read("Structured_Mesh.h5", "/JBV")';
JBP(:, [1,2]) = JBP(:, [1,2]) + 1;
JBV(:, 1) = JBV(:, 1) + 1;

NumPnts_h = size(Points_h, 1);
NumPnts = size(Points, 1);
NumEles = size(Elements_h, 1);

Lx = max(Points(:, 1)) - min(Points(:, 1));
Ly = max(Points(:, 2)) - min(Points(:, 2));

% __________________ Plot mesh __________________
% figure(1)
% subplot(1,2,1)
% patch('vertices', Points_h, 'faces', Elements_h(:, [1,2,3,6,9,8,7,4]), ...
%     'facevertexcdata', zeros(NumPnts_h, 1), 'edgealpha', 1, 'facealpha', 0)
% pbaspect([Lx, Ly, 1])
% subplot(1,2,2)
% patch('vertices', Points, 'faces', Elements, 'facevertexcdata', ...
%     zeros(NumPnts, 1), 'edgealpha', 1, 'facealpha', 0)
% pbaspect([Lx, Ly, 1])
% hold on
% scatter(Points_h(JBV(:, 1), 1), Points_h(JBV(:, 1), 2), 'o')
% -------------------------------------------
syms xi_e real
syms eta_e real

Psi = [1 / 4.0 * (1 - xi_e) * (1 - eta_e), ...
    1 / 4.0 * (1 + xi_e) * (1 - eta_e), ...
    1 / 4.0 * (1 + xi_e) * (1 + eta_e), ...
    1 / 4.0 * (1 - xi_e) * (1 + eta_e)]';
Phi = [1 / 4. * xi_e * eta_e * (xi_e - 1) * (eta_e - 1), ...
    1 / 2. * eta_e * (1 - xi_e^2) * (eta_e - 1), ...
    1 / 4. * xi_e * eta_e * (xi_e + 1) * (eta_e - 1), ...
    1 / 2. * xi_e * (1 - eta_e^2) * (xi_e - 1), ...
    (1 - xi_e^2) * (1 - eta_e^2), ...
    1 / 2. * xi_e * (1 - eta_e^2) * (xi_e + 1), ...
    1 / 4. * xi_e * eta_e * (xi_e - 1) * (eta_e + 1), ...
    1 / 2. * eta_e * (1 - xi_e^2) * (eta_e + 1), ...
    1 / 4. * xi_e * eta_e * (xi_e + 1) * (eta_e + 1)]';

Phi_xi = diff(Phi, xi_e);
Phi_eta = diff(Phi, eta_e);

xi = [0.932469514203152, 0.661209386466265, 0.238619186083197, -0.932469514203152, -0.661209386466265, -0.238619186083197];
w = [0.171324492379170, 0.360761573048139, 0.467913934572691, 0.171324492379170, 0.360761573048139, 0.467913934572691];

m = 1e3;
Step0;

for n = 0.9:-0.2:0.5
    steps = 0;
    R_absolute = 1e4;
    R_relative = 1e4;

    while steps < 20 && R_absolute > 1e-3 && R_relative > 1e-3
        StepMu;

        disp(['R_ab = ', num2str(R_absolute), ', R_re = ', num2str(R_relative), ', step = ', num2str(steps), ', n = ', num2str(n)])
        steps = steps + 1;

    end

    figure(2)
    subplot(1, 2, 1)
    title("Pressure")
    patch('vertices', Points, 'faces', Elements, 'facevertexcdata', pressure, ...
        'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1);
    hold on
    colorbar;
    pbaspect([Lx, Ly, 1])
    
    figure(2)
    subplot(1, 2, 2)
    title("Velocity vectors")
    patch('vertices', Points_h, 'faces', Elements_h(:, [1, 2, 3, 6, 9, 8, 7, 4]), ...
        'facevertexcdata', zeros(NumPnts_h, 1), 'edgealpha', 1, 'facealpha', 0);
    hold on
    quiver(Points_h(:, 1), Points_h(:, 2), u, v, 1);
    hold on
    pbaspect([Lx, Ly, 1])

    pause
    close(2)
end

