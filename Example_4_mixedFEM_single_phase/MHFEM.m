clc
clear all
close all
currentPath = fileparts(mfilename('fullpath'));

load([currentPath, '/Mesh_and_Boundary_Condition.mat']);

M = [
        2., 0., 1., 0., 1., 0.;
        0., 2., 0., 1., 0., 1.;
        1., 0., 2., 0., 1., 0.;
        0., 1., 0., 2., 0., 1.;
        1., 0., 1., 0., 2., 0.;
        0., 1., 0., 1., 0., 2.    
];

NumSepEdge = NumElements * 3;

K = sparse(NumSepEdge + NumElements + NumGlobalEdges, NumSepEdge + NumElements + NumGlobalEdges);
b = sparse(NumSepEdge + NumElements + NumGlobalEdges, 1);

NonFluxEdgeGlobalID = []; 

for i = 1:NumElements
    ID_1 = Elements(i, 1);
    ID_2 = Elements(i, 2);
    ID_3 = Elements(i, 3);

    P1 = Points(ID_1, :)';
    P2 = Points(ID_2, :)';
    P3 = Points(ID_3, :)';

    N = [zeros(2, 1), P1-P2, P1-P3;
         P2-P1, zeros(2, 1), P2-P3;
         P3-P1, P3-P2, zeros(2, 1)];
    A = triangle_area(P1, P2, P3);
    C_T = zeros(3, 3);
    C_T(1, 1) = norm(P3 - P2);
    C_T(2, 2) = norm(P3 - P1);
    C_T(3, 3) = norm(P1 - P2);
    
    B = 1 / 48 / A * C_T' * N' * M * N * C_T;
    
    LID = [ID_2, ID_3, ID_1];
    for j = 1:3

        ID_L1 = LID(j);
        ID_L2 = LID(mod(j, 3) + 1);
        
        Length_pp = norm(Points(ID_L1, :) - Points(ID_L2, :));

        if AttributeEdge(i, mod(j, 3) + 1) == 2 % internal
            GlobalEdgeIDinternal = InternalEdgeNumbering(ID_L1, ID_L2);
            K((i-1) * 3 + j, NumSepEdge + NumElements + GlobalEdgeIDinternal) = -Length_pp;
            K(NumSepEdge + NumElements + GlobalEdgeIDinternal, (i-1) * 3 + j) = -Length_pp;
        end

        if AttributeEdge(i, mod(j, 3) + 1) == -1 % non-flux
            %NonFluxEdgeGlobalID = [NonFluxEdgeGlobalID; GlobalEdgeID];
            B(j, :) = 0;
            B(:, j) = 0;
            B(j, j) = 1;
            C_T(j, j) = 0;
        end

        if AttributeEdge(i, mod(j, 3) + 1) == -2 % outlet
            b((i - 1) * 3 + j) = b((i - 1) * 3 + j) + 1 * Length_pp;
        end

        if AttributeEdge(i, mod(j, 3) + 1) == -3 % outlet
            b((i - 1) * 3 + j) = b((i - 1) * 3 + j) + 100 * Length_pp;
        end
    end

    K((i-1) * 3 + 1:i*3, (i-1) * 3 + 1:i*3) = B;
    K((i-1) * 3 + 1:i*3, NumSepEdge + i) = diag(C_T);
    K(NumSepEdge + i, (i-1) * 3 + 1:i*3) = diag(C_T)';
end

% K = triu(K) + triu(K, 1).'; % Copy the upper triangle to the lower triangle

% for i = size(NonFluxEdgeGlobalID, 1) % non-flux boundary
%     NonFluxEdgeGlobalID_local = NonFluxEdgeGlobalID(i);
% 
%     b = b - K(:, 2 * NumSepEdge + NonFluxEdgeGlobalID_local) .* 0;
% 
%     K(2 * NumSepEdge + NonFluxEdgeGlobalID_local, :) = 0;
%     K(:, 2 * NumSepEdge + NonFluxEdgeGlobalID_local) = 0;
% 
%     K(2 * NumSepEdge + NonFluxEdgeGlobalID_local, 2 * NumSepEdge + NonFluxEdgeGlobalID_local) = 1;
% end

x = inv(K) *  b;

ElementPressue = full(x(NumSepEdge + 1: NumSepEdge + NumElements));
NormalvelocityEdge = full(x(1:NumSepEdge));

figure(1); title('Pressure', 'interpreter', 'latex'); hold on
patch('Vertices', Points, 'Faces', Elements, 'FaceVertexCData', ElementPressue, 'FaceColor', 'flat', 'EdgeAlpha', 0.9, 'facealpha', 1);
hold on
pbaspect([1, 1, 1])
colorbar
caxis([1, 100])


% CentersEdge = zeros(NumSepEdge, 2);
% VeloEdge = zeros(NumSepEdge, 2);
% for i = 1:NumElements
%     ID_1 = Elements(i, 1);
%     ID_2 = Elements(i, 2);
%     ID_3 = Elements(i, 3);
% 
%     P1 = Points(ID_1, :);
%     P2 = Points(ID_2, :);
%     P3 = Points(ID_3, :);
% 
%     Coord = [P2; P3; P1];
% 
%     for j = 1:3
%         V = Coord(mod(j, 3) + 1, :) - Coord(j, :);
%         V = [V(2), -V(1)];
% 
%         V = V ./ norm(V);
% 
%         V = V .* NormalvelocityEdge((i - 1) * 3 + j);
%         VeloEdge((i - 1) * 3 + j, :) = -V;
%         CentersEdge((i - 1) * 3 + j, :) = 0.5*(Coord(mod(j, 3) + 1, :) + Coord(j, :));
%     end
% end
% quiver(CentersEdge(:, 1), CentersEdge(:, 2), VeloEdge(:, 1), VeloEdge(:, 2), 'color', 'r'); hold on

q_ele = zeros(2, NumElements);
center_grid = zeros(NumElements, 2);

InFlux = 0;
OutFlux = 0;

for i = 1:NumElements

    ID_1 = Elements(i, 1);
    ID_2 = Elements(i, 2);
    ID_3 = Elements(i, 3);

    P1 = Points(ID_1, :)';
    P2 = Points(ID_2, :)';
    P3 = Points(ID_3, :)';

    center_s = [(P1 + P2 + P3) ./ 3];
    center_grid(i, :) = center_s';
    A = triangle_area(P1, P2, P3);

    LID = [ID_2, ID_3, ID_1];
    for j = 1:3
        ID_L1 = LID(j);
        ID_L2 = LID(mod(j, 3) + 1);
        Length_pp = norm(Points(ID_L1, :) - Points(ID_L2, :));
        q_ele(:, i) = q_ele(:, i) - Length_pp / A * (center_s - Points(Elements(i, j), :)') .* NormalvelocityEdge((i - 1) * 3 + j);

        if (AttributeEdge(i, mod(j, 3) + 1) == -3) % inlet
            InFlux = InFlux - Length_pp * NormalvelocityEdge((i - 1) * 3 + j);
        elseif (AttributeEdge(i, mod(j, 3) + 1) == -2) % outlet
            OutFlux = OutFlux + Length_pp * NormalvelocityEdge((i - 1) * 3 + j);
        end
    end
end
q_ele = q_ele';
hold on
quiver(center_grid(:, 1), center_grid(:, 2), q_ele(:, 1), q_ele(:, 2), 'color', 'r'); hold on
InFlux
OutFlux