clear;

% %coordinates of the triangle vertex (x1, y1, x2, y2, etc) (random)
% coordElem = [ 2 2; 7 1; 5 4];

%coordinates of the quadrilater vertex (x1, y1, x2, y2, etc) (random)
coordElem = [ 1 1; 3 2; 5 4; 2 5];


if length(coordElem) == 3

    J = triangleElement(coordElem).assembleJ();

elseif length(coordElem) == 4

    J = quadrilaterElement(coordElem).assembleJ();

else

    disp("Wrong coordinates matrix");

end

detJ = det(J);