mu = m;

Dims = 2 * NumPnts_h + NumPnts;

K = sparse(Dims, Dims);
b = sparse(Dims, 1);

pnt_x_9 = reshape(Points_h(Elements_h(:), 1), [16, 9])';
pnt_y_9 = reshape(Points_h(Elements_h(:), 2), [16, 9])';
pnt_x_4 = reshape(Points_h(Elements(:), 1), [16, 4])';
pnt_y_4 = reshape(Points_h(Elements(:), 2), 16, 4)';

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

