clc
clear all
close all

Points = h5read("Structured_Mesh.h5", "/Points")';
Points_h = h5read("Structured_Mesh.h5", "/Points_h")';
Elements = h5read("Structured_Mesh.h5", "/Elements")' + 1;
Elements_h = h5read("Structured_Mesh.h5", "/Elements_h")' + 1;

NumPnts_h = size(Points_h, 1);
NumPnts = size(Points, 1);
NumEles = size(Elements_h, 1);

Lx = max(Points(:, 1)) - min(Points(:, 1));
Ly = max(Points(:, 2)) - min(Points(:, 2));
% —————————— re-number the node______________
Visited = zeros(NumPnts_h, 1);
Points_new = zeros(NumPnts_h, 2);
Points_new(1:NumPnts, :) = Points(:, 1:2);

TmpPnt = NumPnts + 1;
Elements_new = Elements_h(:, [2, 4, 5, 6, 8]);

for i = 1:NumEles
    for j = 1:5
        Pnt_ID_y = Elements_new(i, j);
        if Visited(Pnt_ID_y) == 0
            Visited(Pnt_ID_y) = TmpPnt;
            Points_new(TmpPnt, :) = Points_h(Elements_new(i, j), [1, 2]);
            Elements_new(i, j) = TmpPnt;
            TmpPnt = TmpPnt + 1;
        else
            Elements_new(i, j) = Visited(Pnt_ID_y);
        end
    end
end
Elements_h(:, [1, 3, 9, 7]) = Elements;
Elements_h(:, [2, 4, 5, 6, 8]) = Elements_new;
Points_h = Points_new;
clear Elements_new Points_new TmpPnt Pnt_ID_y Visited i j

% __________________ Plot mesh __________________
figure(1)
subplot(1,2,1)
patch('vertices', Points_h, 'faces', Elements_h(:, [1,2,3,6,9,8,7,4]), ...
    'facevertexcdata', zeros(NumPnts_h, 1), 'edgealpha', 1, 'facealpha', 0)
pbaspect([Lx, Ly, 1])
subplot(1,2,2)
patch('vertices', Points, 'faces', Elements, 'facevertexcdata', ...
    zeros(NumPnts, 1), 'edgealpha', 1, 'facealpha', 0)
pbaspect([Lx, Ly, 1])


% _________________ solver _______________
mu = 1e3;
xi = [0.932469514203152, 0.661209386466265, 0.238619186083197, -0.932469514203152, -0.661209386466265, -0.238619186083197];
w = [0.171324492379170, 0.360761573048139, 0.467913934572691, 0.171324492379170, 0.360761573048139, 0.467913934572691];

Dims = 2 * NumPnts_h + NumPnts;

K = sparse(Dims, Dims);
b = sparse(Dims, 1);

velocity_condition = zeros(NumPnts_h, 3); % ID, u, v
velocity_condition_counter = 1;

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
for e = 1:NumEles
    disp([num2str(e), '/', num2str(NumEles)])
    D11 = zeros(9, 9);
    D12 = zeros(9, 9);
    D21 = zeros(9, 9);
    D22 = zeros(9, 9);

    C1 = zeros(9, 4);
    C2 = zeros(9, 4);

    B1 = zeros(4, 9);
    B2 = zeros(4, 9);

    F1 = zeros(9, 1);
    F2 = zeros(9, 1);

    pnt_x_9 = Points_h(Elements_h(e, :), 1);
    pnt_y_9 = Points_h(Elements_h(e, :), 2);
    pnt_x_4 = Points_h(Elements(e, :), 1);
    pnt_y_4 = Points_h(Elements(e, :), 2);

    for i = 1:6
        xi_e_i = xi(i);

        for j = 1:6

            eta_e_j = xi(j);

            Phi_xi_e = double(subs(Phi_xi, [xi_e, eta_e], [xi_e_i, eta_e_j]));
            Phi_eta_e = double(subs(Phi_eta, [xi_e, eta_e], [xi_e_i, eta_e_j]));

            x_xi = dot(Phi_xi_e, pnt_x_9);
            y_xi = dot(Phi_xi_e, pnt_y_9);
            x_eta = dot(Phi_eta_e, pnt_x_9);
            y_eta = dot(Phi_eta_e, pnt_y_9);

            J = [x_xi, y_xi; x_eta, y_eta];
            AAAA = J \ [Phi_xi_e'; Phi_eta_e'];
            Phi_x = AAAA(1, :)';
            Phi_y = AAAA(2, :)';
            detJ = det(J);

            D11 = D11 + mu * w(i) * w(j) * (2 * (Phi_x * Phi_x') + (Phi_y * Phi_y')) * detJ;
            D12 = D12 + mu * w(i) * w(j) * ((Phi_x * Phi_y')) * detJ;
            D21 = D21 + mu * w(i) * w(j) * ((Phi_y * Phi_x')) * detJ;
            D22 = D22 + mu * w(i) * w(j) * (2 * (Phi_y * Phi_y') + (Phi_x * Phi_x')) * detJ;

            Psi_e = double(subs(Psi, [xi_e, eta_e], [xi_e_i, eta_e_j]));

            C1 = C1 + w(i) * w(j) * (Phi_x * Psi_e') * detJ;
            C2 = C2 + w(i) * w(j) * (Phi_y * Psi_e') * detJ;

            B1 = B1 + w(i) * w(j) * (Psi_e * Phi_x') * detJ;
            B2 = B2 + w(i) * w(j) * (Psi_e * Phi_y') * detJ;
        end
    end

    for l = 1:4 %%%% boundary condition
        if pnt_x_4(l) == 0 && pnt_x_4(mod(l, 4)+1) == 0 || pnt_x_4(l) == 80 && pnt_x_4(mod(l, 4)+1) == 80
            p_i_e = zeros(4, 1);
            Length_boundary = sqrt((pnt_x_4(l) - pnt_x_4(mod(l, 4)+1))^2+ ...
                (pnt_y_4(l) - pnt_y_4(mod(l, 4)+1))^2);
            if pnt_x_4(l) == 0 && pnt_x_4(mod(l, 4)+1) == 0
                p_i_e([l, mod(l, 4) + 1]) = -1000;
            elseif pnt_x_4(l) == .80 && pnt_x_4(mod(l, 4)+1) == .80
                p_i_e([l, mod(l, 4) + 1]) = 0;
            end
            cos_theta_x = 1;
            cos_theta_y = 0;
            for i = 1:6
                xi_e_i = xi(i);
                eta_e_i = xi(i);
                if l == 1
                    eta_e_i = -1;
                elseif l == 2
                    xi_e_i = 1;
                elseif l == 3
                    eta_e_i = 1;
                elseif l == 4
                    xi_e_i = -1;
                end
                Phi_e = double(subs(Phi, [xi_e, eta_e], [xi_e_i, eta_e_i]));
                Psi_e = double(subs(Psi, [xi_e, eta_e], [xi_e_i, eta_e_i]));
                F1 = F1 + w(i) * Phi_e * Psi_e' * p_i_e * cos_theta_x * Length_boundary / 2;
                F2 = F2 + w(i) * Phi_e * Psi_e' * p_i_e * cos_theta_y * Length_boundary / 2;
            end
        elseif pnt_y_4(l) == 0 && pnt_y_4(mod(l, 4)+1) == 0 ...
                || pnt_y_4(l) == .6 && pnt_y_4(mod(l, 4)+1) == .6

            ID_reOrgnized = Elements_h(e, [1, 2, 3, 6, 9, 8, 7, 4]);
            localID_v = (l - 1) * 2 + 1;
            %%%[pnt_y_4(l), pnt_y_4(mod(l, 4)+1)]
            if pnt_y_4(l) == 0 && pnt_y_4(mod(l, 4)+1) == 0

                velocity_condition(velocity_condition_counter, :) = [double(ID_reOrgnized(localID_v)), 0, 0];
                velocity_condition(velocity_condition_counter+1, :) = [double(ID_reOrgnized(localID_v+1)), 0, 0];
                velocity_condition(velocity_condition_counter+2, :) = [double(ID_reOrgnized(mod(localID_v+2, 8))), 0, 0];
                velocity_condition_counter = velocity_condition_counter + 3;
            else
                velocity_condition(velocity_condition_counter, :) = [double(ID_reOrgnized(localID_v)), 1e-2, 0];
                velocity_condition(velocity_condition_counter+1, :) = [double(ID_reOrgnized(localID_v+1)), 1e-2, 0];
                velocity_condition(velocity_condition_counter+2, :) = [double(ID_reOrgnized(mod(localID_v+2, 8))), 1e-2, 0];
                velocity_condition_counter = velocity_condition_counter + 3;
            end
        end
    end

    % ______________ assemble __________
    K(Elements_h(e, :), Elements_h(e, :)) = K(Elements_h(e, :), Elements_h(e, :)) + D11;
    K(Elements_h(e, :), Elements_h(e, :)+NumPnts_h) = K(Elements_h(e, :), Elements_h(e, :)+NumPnts_h) + D12;
    K(Elements_h(e, :)+NumPnts_h, Elements_h(e, :)) = K(Elements_h(e, :)+NumPnts_h, Elements_h(e, :)) + D21;
    K(Elements_h(e, :)+NumPnts_h, Elements_h(e, :)+NumPnts_h) = K(Elements_h(e, :)+NumPnts_h, Elements_h(e, :)+NumPnts_h) + D22;

    K(Elements_h(e, :), Elements(e, :)+2*NumPnts_h) = K(Elements_h(e, :), Elements(e, :)+2*NumPnts_h) - C1;
    K(Elements_h(e, :)+NumPnts_h, Elements(e, :)+2*NumPnts_h) = ...
        K(Elements_h(e, :)+NumPnts_h, Elements(e, :)+2*NumPnts_h) - C2;

    K(Elements(e, :)+2*NumPnts_h, Elements_h(e, :)) = K(Elements(e, :)+2*NumPnts_h, Elements_h(e, :)) + ...
        B1;
    K(Elements(e, :)+2*NumPnts_h, Elements_h(e, :)+NumPnts_h) = K(Elements(e, :)+2*NumPnts_h, Elements_h(e, :)+NumPnts_h) + ...
        B2;

    b(Elements_h(e, :), 1) = b(Elements_h(e, :), 1) - F1;
    b(Elements_h(e, :)+NumPnts_h, 1) = b(Elements_h(e, :)+NumPnts_h, 1) - F2;
end

velocity_condition(find(velocity_condition(:, 1) == 0, 1):end, :) = [];

% ___________ velocity boundary _______________
for i = 1:size(velocity_condition, 1)
    PntID_v = velocity_condition(i, 1);

    b = b - K(:, PntID_v) .* velocity_condition(i, 2);
    K(:, PntID_v) = 0;
    K(PntID_v, :) = 0;
    K(PntID_v, PntID_v) = 1;
    b(PntID_v, 1) = velocity_condition(i, 2);

    PntID_v = PntID_v + NumPnts_h;
    b = b - K(:, PntID_v) .* velocity_condition(i, 3);
    K(:, PntID_v) = 0;
    K(PntID_v, :) = 0;
    K(PntID_v, PntID_v) = 1;
    b(PntID_v, 1) = velocity_condition(i, 3);
end

x = K \ b;

pressure = full(x(2*NumPnts_h+1:end));
u = full(x(1:NumPnts_h));
v = full(x(NumPnts_h+1:2*NumPnts_h));

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