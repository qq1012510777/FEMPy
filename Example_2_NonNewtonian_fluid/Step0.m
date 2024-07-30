mu = m;

Dims = 2 * NumPnts_h + NumPnts;

K = sparse(Dims, Dims);
b = sparse(Dims, 1);

pnt_x_9 = reshape(Points_h(Elements_h(:), 1), [NumEles, 9])';
pnt_y_9 = reshape(Points_h(Elements_h(:), 2), [NumEles, 9])';
pnt_x_4 = reshape(Points_h(Elements(:), 1), [NumEles, 4])';
pnt_y_4 = reshape(Points_h(Elements(:), 2), NumEles, 4)';

D11 = zeros(9, 9, NumEles);
D12 = zeros(9, 9, NumEles);
D21 = zeros(9, 9, NumEles);
D22 = zeros(9, 9, NumEles);

C1 = zeros(9, 4, NumEles);
C2 = zeros(9, 4, NumEles);

B1 = zeros(4, 9, NumEles);
B2 = zeros(4, 9, NumEles);

for i = 1:6
    xi_e_i = xi(i);
    for j = 1:6
        eta_e_j = xi(j);

        Phi_xi_e = double(subs(Phi_xi, [xi_e, eta_e], [xi_e_i, eta_e_j]));
        Phi_eta_e = double(subs(Phi_eta, [xi_e, eta_e], [xi_e_i, eta_e_j]));
        Psi_e = double(subs(Psi, [xi_e, eta_e], [xi_e_i, eta_e_j]));

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

        D11 = D11 + mu * w(i) * w(j) * pagemtimes((2 * pagemtimes(Phi_x, reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)])) ...
            +pagemtimes(Phi_y, reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)]))), DetJ);
        D12 = D12 + mu * w(i) * w(j) * pagemtimes((pagemtimes(Phi_x, reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)])) ...
            ), DetJ);
        D21 = D21 + mu * w(i) * w(j) * pagemtimes((pagemtimes(Phi_y, reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)])) ...
            ), DetJ);
        D22 = D22 + mu * w(i) * w(j) * pagemtimes((2 * pagemtimes(Phi_y, reshape(Phi_y, [size(Phi_y, 2), size(Phi_y, 1), size(Phi_y, 3)])) ...
            +pagemtimes(Phi_x, reshape(Phi_x, [size(Phi_x, 2), size(Phi_x, 1), size(Phi_x, 3)]))), DetJ);

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
    B1(:); ...
    B2(:)]];

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

        if JBP(e, 2) == 1
            eta_e_i = -1;
        elseif JBP(e, 2) == 2
            xi_e_i = 1;
        elseif JBP(e, 2) == 3
            eta_e_i = 1;
        elseif JBP(e, 2) == 4
            xi_e_i = -1;
        end

        Phi_e = double(subs(Phi, [xi_e, eta_e], [xi_e_i, eta_e_i]));
        Psi_e = double(subs(Psi, [xi_e, eta_e], [xi_e_i, eta_e_i]));

        F1 = F1 + sign(cos_theta_x) * w(i) * Phi_e * Psi_e' * p_i_e * cos_theta_x * Length_boundary / 2;
        F2 = F2 + sign(cos_theta_y) * w(i) * Phi_e * Psi_e' * p_i_e * cos_theta_y * Length_boundary / 2;
    end

    b(Elements_h(eleID, :), 1) = b(Elements_h(e, :), 1) - F1;
    b(Elements_h(eleID, :) + NumPnts_h, 1) = b(Elements_h(e, :)+NumPnts_h, 1) - F2;

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