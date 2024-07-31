Dims = 2 * NumPnts_h + NumPnts;

K = sparse(Dims, Dims);
b = sparse(Dims, 1);

pnt_x_9 = reshape(Points_h(Elements_h(:), 1), [NumEles, 9])';
pnt_y_9 = reshape(Points_h(Elements_h(:), 2), [NumEles, 9])';
pnt_x_4 = reshape(Points(Elements(:), 1), [NumEles, 4])';
pnt_y_4 = reshape(Points(Elements(:), 2), [NumEles, 4])';

D11 = zeros(9, 9, NumEles);
D12 = zeros(9, 9, NumEles);
D21 = zeros(9, 9, NumEles);
D22 = zeros(9, 9, NumEles);

C1 = zeros(9, 4, NumEles);
C2 = zeros(9, 4, NumEles);

B1 = zeros(4, 9, NumEles);
B2 = zeros(4, 9, NumEles);

Mu_ElementWise = reshape(Mu_eachPnt(Elements_h'), 9, 1, NumEles);

for i = 1:6
    xi_e_i = xi(i);
    for j = 1:6
        eta_e_j = xi(j);

        Phi_xi_e = double(subs(Phi_xi, [xi_e, eta_e], [xi_e_i, eta_e_j]));
        Phi_eta_e = double(subs(Phi_eta, [xi_e, eta_e], [xi_e_i, eta_e_j]));
        Psi_e = double(subs(Psi, [xi_e, eta_e], [xi_e_i, eta_e_j]));
        Phi_e = double(subs(Phi, [xi_e, eta_e], [xi_e_i, eta_e_j]));

        Mu_h = pagemtimes(Phi_e', Mu_ElementWise);

        x_xi = Phi_xi_e' * pnt_x_9;
        y_xi = Phi_xi_e' * pnt_y_9;
        x_eta = Phi_eta_e' * pnt_x_9;
        y_eta = Phi_eta_e' * pnt_y_9;

        J = zeros(2, 2, NumEles);
        J(1, 1, :) = x_xi;
        J(1, 2, :) = y_xi;
        J(2, 1, :) = x_eta;
        J(2, 2, :) = y_eta;

        AAAA = pagemldivide(J, [Phi_xi_e'; Phi_eta_e']);

        Phi_x = AAAA(1, :, :);
        Phi_y = AAAA(2, :, :);

        Phi_x = reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)]);
        Phi_y = reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)]);

        EigenvalueJ = pageeig(J);
        DetJ = EigenvalueJ(1, 1, :) .* EigenvalueJ(2, 1, :);

        % D11 = D11 + mu * w(i) * w(j) * (2 * (Phi_x * Phi_x') + (Phi_y * Phi_y')) * detJ;

        D11 = D11 + w(i) * w(j) * pagemtimes(Mu_h, pagemtimes((2 * pagemtimes(Phi_x, reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)])) ...
            +pagemtimes(Phi_y, reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)]))), DetJ));
        D12 = D12 + w(i) * w(j) * pagemtimes(Mu_h, pagemtimes((pagemtimes(Phi_x, reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)])) ...
            ), DetJ));
        D21 = D21 + w(i) * w(j) * pagemtimes(Mu_h, pagemtimes((pagemtimes(Phi_y, reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)])) ...
            ), DetJ));
        D22 = D22 + w(i) * w(j) * pagemtimes(Mu_h, pagemtimes((2 * pagemtimes(Phi_y, reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)])) ...
            +pagemtimes(Phi_x, reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)]))), DetJ));

        % C1 = C1 + w(i) * w(j) * (Phi_x * Psi_e') * detJ;
        C1 = C1 + w(i) * w(j) * pagemtimes(pagemtimes(Phi_x, Psi_e'), DetJ);
        C2 = C2 + w(i) * w(j) * pagemtimes(pagemtimes(Phi_y, Psi_e'), DetJ);

        % B1 = B1 + w(i) * w(j) * (Psi_e * Phi_x') * detJ;
        B1 = B1 + w(i) * w(j) * pagemtimes(pagemtimes(Psi_e, reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)])), DetJ);
        B2 = B2 + w(i) * w(j) * pagemtimes(pagemtimes(Psi_e, reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)])), DetJ);

    end
end

% _____________ D index ____________
Idx_h_row = reshape(Elements_h', [9, 1, NumEles]);
Idx_h_row = repmat(Idx_h_row, 1, 9, 1);
Idx_h_col = reshape(Elements_h', [9, 1, NumEles]);
Idx_h_col = pagetranspose(Idx_h_col);
Idx_h_col = repmat(Idx_h_col, 9, 1, 1);

D11_Idx = (Idx_h_row - 1) * Dims + Idx_h_col;
D12_Idx = (Idx_h_row - 1) * Dims + Idx_h_col + NumPnts_h;
D21_Idx = (Idx_h_row - 1 + NumPnts_h) * Dims + Idx_h_col;
D22_Idx = (Idx_h_row - 1 + NumPnts_h) * Dims + Idx_h_col + NumPnts_h;
%
% _____________ C index ____________
Idx_h_row = reshape(Elements_h', [9, 1, NumEles]);
Idx_h_row = repmat(Idx_h_row, 1, 4, 1);
Idx_h_col = reshape(Elements', [4, 1, NumEles]);
Idx_h_col = pagetranspose(Idx_h_col);
Idx_h_col = repmat(Idx_h_col, 9, 1, 1);

C1_Idx = (Idx_h_row - 1) * Dims + Idx_h_col + 2 * NumPnts_h;
C2_Idx = (Idx_h_row - 1 + NumPnts_h) * Dims + Idx_h_col + 2 * NumPnts_h;

% _____________ B index ____________
Idx_h_row = reshape(Elements', [4, 1, NumEles]);
Idx_h_row = repmat(Idx_h_row, 1, 9, 1);
Idx_h_col = reshape(Elements_h', [9, 1, NumEles]);
Idx_h_col = pagetranspose(Idx_h_col);
Idx_h_col = repmat(Idx_h_col, 4, 1, 1);

B1_Idx = (Idx_h_row - 1 + 2 * NumPnts_h) * Dims + Idx_h_col;
B2_Idx = (Idx_h_row - 1 + 2 * NumPnts_h) * Dims + Idx_h_col + NumPnts_h;

IndexValue = [double([D11_Idx(:); ...
    D12_Idx(:); ...
    D21_Idx(:); ...
    D22_Idx(:); ...
    C1_Idx(:); ...
    C2_Idx(:); ...
    B1_Idx(:); ...
    B2_Idx(:)]), [D11(:); ...
    D12(:); ...
    D21(:); ...
    D22(:); ...
    -C1(:); ...
    -C2(:); ...
    -B1(:); ...
    -B2(:)]];

IndexValue = sortrows(IndexValue, 1);

% sort the elements in the same positions, sum them if they are in the same positions
[A2u, ~, ix] = unique(IndexValue(:, 1));
A1sums = accumarray(ix, IndexValue(:, 2), [], @sum);
IndexValue = [A2u, A1sums];

K(IndexValue(:, 1)) = IndexValue(:, 2);

% _________________ boundary condition Pressure

for e = 1:size(JBP, 1)
    p_i_e = zeros(4, 1);

    l = JBP(e, 2);
    eleID = JBP(e,1);
    pnt_x_4 = Points_h(Elements(eleID, :), 1);
    pnt_y_4 = Points_h(Elements(eleID, :), 2); 

    Length_boundary = sqrt((pnt_x_4(l) - pnt_x_4(mod(l, 4)+1))^2+ ...
        (pnt_y_4(l) - pnt_y_4(mod(l, 4)+1))^2);
    p_i_e([l, mod(l, 4) + 1]) = [JBP(e,3), JBP(e,4)];

    cos_theta_x = JBP(e,5);
    cos_theta_y = JBP(e,6);
    
    F1 = zeros(9, 1);
    F2 = zeros(9, 1);

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

    b(Elements_h(eleID, :), 1) = b(Elements_h(e, :), 1) - F1;
    b(Elements_h(eleID, :) + NumPnts_h, 1) = b(Elements_h(e, :) + NumPnts_h, 1) - F2;

end

% ____________boundary condition velocity
for i = 1:size(JBV, 1)
    PntID_v = JBV(i, 1);
    
    b = b - K(:, PntID_v) .* JBV(i, 2);
    K(:, PntID_v) = 0;
    K(PntID_v, :) = 0;
    K(PntID_v, PntID_v) = 1;
    b(PntID_v, 1) = JBV(i, 2);

    PntID_v = PntID_v + NumPnts_h;
    b = b - K(:, PntID_v) .* JBV(i, 3);
    K(:, PntID_v) = 0;
    K(PntID_v, :) = 0;
    K(PntID_v, PntID_v) = 1;
    b(PntID_v, 1) = JBV(i, 3);
end

x = K \ b;
pressure_1 = full(x(2*NumPnts_h+1:end));
u_1 = full(x(1:NumPnts_h));
v_1 = full(x(NumPnts_h+1:2*NumPnts_h));

clear K b A1sums A2u AAAA D22 D11 D12 D21 
clear D11_Idx D12_Idx D21_Idx  D22_Idx   C1_Idx   C2_Idx    B1_Idx   B2_Idx
clear B1 B2 C1 C2 cos_theta_y cos_theta_x DetJ Dims e EigenvalueJ
clear eleID eta_e_i eta_e_j F1 F2 i Idx_h_col Idx_h_row IndexValue
clear ix j J eta l  Length_boundary p_i_e Phi_e Phi_eta_e Phi_x 
clear Phi_xi_e Phi_y pnt_x_4 pnt_y_4 PntID_v
clear Psi_e x_eta x_xi xi_ei y_eta y_xi xi_e_i mu


% figure(2)
% subplot(1, 2, 1)
% title("Pressure")
% patch('vertices', Points, 'faces', Elements, 'facevertexcdata', pressure_1, ...
%     'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1);
% hold on
% colorbar;
% pbaspect([Lx, Ly, 1])
% 
% figure(2)
% subplot(1, 2, 2)
% title("Velocity vectors")
% patch('vertices', Points_h, 'faces', Elements_h(:, [1, 2, 3, 6, 9, 8, 7, 4]), ...
%     'facevertexcdata', zeros(NumPnts_h, 1), 'edgealpha', 1, 'facealpha', 0);
% hold on
% quiver(Points_h(:, 1), Points_h(:, 2), u_1, v_1, 1);
% hold on
% pbaspect([Lx, Ly, 1])


% ____________________________calculate mu (non-Newtonian) ___________-
ElementWise_u = reshape(u(Elements_h(:), 1), [9, 1, NumEles]);
ElementWise_v = reshape(v(Elements_h(:), 1), [9, 1, NumEles]);

xi_mu = [-1, 0, 1, -1, 0, 1, -1, 0, 1]';
eta_mu = [-1, -1, -1, 0, 0, 0, 1, 1, 1]';

mu_dynamic = zeros(NumEles, 9);
AreaElements = zeros(NumEles, 9);
for i = 1:9
    Phi_xi_e = double(subs(Phi_xi, [xi_e, eta_e], [xi_mu(i), eta_mu(i)]));
    Phi_eta_e = double(subs(Phi_eta, [xi_e, eta_e], [xi_mu(i), eta_mu(i)]));

    x_xi = Phi_xi_e' * pnt_x_9;
    y_xi = Phi_xi_e' * pnt_y_9;
    x_eta = Phi_eta_e' * pnt_x_9;
    y_eta = Phi_eta_e' * pnt_y_9;

    J = zeros(2, 2, NumEles);
    J(1, 1, :) = x_xi;
    J(1, 2, :) = y_xi;
    J(2, 1, :) = x_eta;
    J(2, 2, :) = y_eta;

    AAAA = pagemldivide(J, [Phi_xi_e'; Phi_eta_e']);

    Phi_x = AAAA(1, :, :);
    Phi_y = AAAA(2, :, :);
    EigenvalueJ = pageeig(J);
    DetJ = EigenvalueJ(1, 1, :) .* EigenvalueJ(2, 1, :);

    AreaElements(:, i) = squeeze(DetJ);

    Phi_x_times_u = squeeze(pagemtimes(Phi_x, ElementWise_u).^2);
    Phi_y_times_v = squeeze(pagemtimes(Phi_y, ElementWise_v).^2);
    Phi_x_times_v = squeeze(pagemtimes(Phi_x, ElementWise_v));
    Phi_y_times_u = squeeze(pagemtimes(Phi_y, ElementWise_u));

    mu_dynamic(:, i) = m .* (2 .* Phi_x_times_u + 2 .* Phi_y_times_v + (Phi_x_times_v + Phi_y_times_u).^2).^((n - 1) / 2.);

end

Mu_eachPnt_1 = [double(Elements_h(:)), mu_dynamic(:) .* AreaElements(:)];
Mu_eachPnt_1 = sortrows(Mu_eachPnt_1, 1);

% sort the elements in the same positions, sum them if they are in the same positions
[A2u, ~, ix] = unique(Mu_eachPnt_1(:, 1));
A1sums = accumarray(ix, Mu_eachPnt_1(:, 2), [], @sum);
Mu_eachPnt_1 = [A2u, A1sums];

Mu_sumArea = [double(Elements_h(:)), AreaElements(:)];
Mu_sumArea = sortrows(Mu_sumArea, 1);

% sort the elements in the same positions, sum them if they are in the same positions
[A2u, ~, ix] = unique(Mu_sumArea(:, 1));
A1sums = accumarray(ix, Mu_sumArea(:, 2), [], @sum);
Mu_sumArea = [A2u, A1sums];

Mu_eachPnt_1(:, 2) = Mu_eachPnt_1(:, 2) ./ Mu_sumArea(:, 2);

Mu_eachPnt_1(:, 1) = [];

Mu_eachPnt(Mu_eachPnt > 1e10) = 1e10;
Mu_eachPnt(Mu_eachPnt < 1e-10) = 1e-10;

clear Mu_sumArea A1sums ix AreaElements ElementWise_u ElementWise_v
clear xi_mu eta_mu mu_dynamic A2u
clear Mu_sumArea A1sums ix AreaElements ElementWise_u ElementWise_v
clear xi_mu eta_mu mu_dynamic A2u
clear K b A1sums A2u AAAA D22 D11 D12 D21 
clear D11_Idx D12_Idx D21_Idx  D22_Idx   C1_Idx   C2_Idx    B1_Idx   B2_Idx
clear B1 B2 C1 C2 cos_theta_y cos_theta_x DetJ Dims e EigenvalueJ
clear eleID eta_e_i eta_e_j F1 F2 i Idx_h_col Idx_h_row IndexValue
clear ix j J eta l  Length_boundary p_i_e Phi_e Phi_eta_e Phi_x 
clear Phi_xi_e Phi_y pnt_x_4 pnt_y_4 PntID_v
clear Psi_e x_eta x_xi xi_ei y_eta y_xi xi_e_i mu
clear AAAA DetJ EigenvalueJ Phi_eta_e Phi_xi_e Phi_x_times_u Phi_y_times_v Phi_x_times_v Phi_y_times_u
clear pnt_x_9 pnt_y_9
%_____________________________________________________________________

R_absolute = max(vecnorm([u', v', pressure', Mu_eachPnt'; u_1', v_1', pressure_1', Mu_eachPnt_1'])) .* 0.01;
R_relative = max(vecnorm([u', v', pressure', Mu_eachPnt'; u_1', v_1', pressure_1', Mu_eachPnt_1']) ...
    ./ abs([u_1', v_1', pressure_1', Mu_eachPnt_1'])) .* 0.01;

u = u_1;
v = v_1;
pressure = pressure_1;
Mu_eachPnt = Mu_eachPnt_1;
% close(2)

clear  u_1 v_1 pressure_1 Mu_eachPnt_1