clear;

%coordinates of the triangle vertex (x1, y1, x2, y2, etc) (random)
coordElem = [ 2 2; 7 1; 5 4];

% %coordinates of the quadrilater vertex (x1, y1, x2, y2, etc) (random)
% coordElem = [ 1 1; 3 2; 5 4; 2 5];

elementShape = length(coordElem);

switch elementShape

    case 3

        [J,N] = TriangleLinearElement(coordElem).assembleJ();

    case 4

        [J,N] = QuadrilaterLinearElement(coordElem).assembleJ();

    otherwise

        disp("Wrong coordinates matrix");

end

detJ = det(J);