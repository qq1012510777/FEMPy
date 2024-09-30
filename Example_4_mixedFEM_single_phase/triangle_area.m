function area = triangle_area(P1, P2, P3)
    % triangleAreaHeron calculates the area of a triangle given three points
    % P1, P2, P3 are 1x2 vectors representing the coordinates of the points
    % Example: P1 = [x1, y1], P2 = [x2, y2], P3 = [x3, y3]

    % Calculate the lengths of the sides
    a = sqrt((P2(1) - P1(1))^2 + (P2(2) - P1(2))^2);
    b = sqrt((P3(1) - P2(1))^2 + (P3(2) - P2(2))^2);
    c = sqrt((P1(1) - P3(1))^2 + (P1(2) - P3(2))^2);

    % Calculate the semi-perimeter
    s = (a + b + c) / 2;

    % Calculate the area using Heron's formula
    area = sqrt(s * (s - a) * (s - b) * (s - c));
end