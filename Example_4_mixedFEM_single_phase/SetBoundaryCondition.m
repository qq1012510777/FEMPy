clc
clear all
close all

currentPath = fileparts(mfilename('fullpath'));

Points = [h5read([currentPath, '/Mesh.h5'], '/Points')]';
Points = Points(:, [1, 2]);
Elements = double([h5read([currentPath, '/Mesh.h5'], '/Elements')]');

NumPoints = size(Points, 1);
NumElements = size(Elements, 1);

MatrixEdge = sparse(NumPoints, NumPoints);

IDtmp = [];
for i = 1:3
    IDtmp = [IDtmp; Elements(:, i) + (Elements(:, mod(i, 3) + 1) - 1) .* NumPoints];
end
MatrixEdge(IDtmp) = 1;

AttributeEdge = Elements .* 0;

%BoundaryEdges = [];
for i = 1:3
    EdgeID = Elements(:, i) + (Elements(:, mod(i, 3) + 1) - 1) .* NumPoints;
    EdgeID_Inv = (Elements(:, i) - 1) .* NumPoints + Elements(:, mod(i, 3) + 1); 
    Ac = full(MatrixEdge(EdgeID) + MatrixEdge(EdgeID_Inv));
    AttributeEdge(:, i) = Ac;  
    %op = find(Ac == 1);
    %BoundaryEdges = [BoundaryEdges; Elements(op, [i, mod(i, 3) + 1])];
end

InternalEdgeNumbering = sparse(NumPoints, NumPoints);

InternalNoEdge = 1;

for i = 1:NumElements
    for j = 1:3
        ID1 = Elements(i, j);
        ID2 = Elements(i, mod(j, 3) + 1);
        
        if (MatrixEdge(ID1, ID2) == 1 && MatrixEdge(ID2, ID1) == 1)
            if (InternalEdgeNumbering(ID1, ID2) == 0 && InternalEdgeNumbering(ID2, ID1) == 0)
                InternalEdgeNumbering(ID1, ID2) = InternalNoEdge;
                InternalNoEdge = InternalNoEdge + 1;
            else
                InternalEdgeNumbering(ID1, ID2) = InternalEdgeNumbering(ID2, ID1);
            end
        end
        
        if AttributeEdge(i, j) ~= 2
            if (Points(ID1, 1) == Points(ID2, 1) && Points(ID1, 1) == 0) % inlet
                AttributeEdge(i, j) = -3;
            elseif (Points(ID1, 1) == Points(ID2, 1) && Points(ID1, 1) == 10) % outlet
                AttributeEdge(i, j) = -2;
            else                                                            % Non-flux
                AttributeEdge(i, j) = -1;
            end
        end
    end
end

BoundaryEdges_inlet = [];
BoundaryEdges_outlet = [];
BoundaryEdges_nonflux = [];

for i = 1:3
    op = find(AttributeEdge(:, i) == -3);
    BoundaryEdges_inlet = [BoundaryEdges_inlet; Elements(op, [i, mod(i, 3) + 1])];

    op = find(AttributeEdge(:, i) == -2);
    BoundaryEdges_outlet = [BoundaryEdges_outlet; Elements(op, [i, mod(i, 3) + 1])];

    op = find(AttributeEdge(:, i) == -1);
    BoundaryEdges_nonflux = [BoundaryEdges_nonflux; Elements(op, [i, mod(i, 3) + 1])];
end

figure(10); title('Check boundary determination', 'interpreter', 'latex'); hold on
patch('Vertices', Points, 'Faces', Elements, 'FaceVertexCData', zeros(size(Points, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 0.9, 'facealpha', 0);
hold on
%patch('Vertices', Points, 'Faces', BoundaryEdges, 'FaceVertexCData', zeros(size(Points, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'r');
patch('Vertices', Points, 'Faces', BoundaryEdges_inlet, 'FaceVertexCData', zeros(size(Points, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'r');hold on
patch('Vertices', Points, 'Faces', BoundaryEdges_outlet, 'FaceVertexCData', zeros(size(Points, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'b'); hold on
patch('Vertices', Points, 'Faces', BoundaryEdges_nonflux, 'FaceVertexCData', zeros(size(Points, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0, 'edgecolor', 'g');hold on
hold on

for i = 1:InternalNoEdge-1
    [row_t, col_t] = find(InternalEdgeNumbering == i);
    KO = 'r';
    if size(row_t, 1) == 1
        KO = 'b';
    end
    
    ID_1 = row_t(1);
    ID_2 = col_t(1);

    Centers = (Points(ID_1, :) + Points(ID_2, :)) .* 0.5;
    text(Centers(1), Centers(2), num2str(i), 'color', KO); hold on
end

hold on
pbaspect([1, 1, 1])
NumGlobalEdges = InternalNoEdge-1;
save([currentPath, '/Mesh_and_Boundary_Condition.mat'], 'Points', 'Elements', 'AttributeEdge', 'InternalEdgeNumbering', 'NumPoints', 'NumElements', 'NumGlobalEdges');