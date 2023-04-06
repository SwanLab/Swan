clear;

%coordinates of the triangle vertex (x1,y1; x2,y2; etc) (random)
coord = [2 2; 7 1; 5 4];

% xi and eta equation, found on article Chapter 4 Finite Element
% Approximation page 6
detDenXi = det([1 coord(2,1) coord(2,2); 1 coord(1,1) coord(1,2); 1 coord(3,1) coord(3,2)]);
detDenEta = det([1 coord(3,1) coord(3,2); 1 coord(1,1) coord(1,2); 1 coord(2,1) coord(2,2)]);

coefXi = [coord(1,2)-coord(3,2), coord(3,1)-coord(1,1),coord(1,1)*coord(3,2)-coord(1,2)*coord(3,1)] / detDenXi;
coefEta = [coord(1,2)-coord(2,2), coord(2,1)-coord(1,1),coord(1,1)*coord(2,2)-coord(1,2)*coord(2,1)] / detDenEta;

coef = [coefXi; coefEta; 1-coefXi-coefEta]; %Matrix of coefficients

% AD using my code.
x = ValDerForward(1,[1 0]);
y = ValDerForward(1,[0 1]);

n1 = coef(1,1)*x + coef(1,2)*y + coef(1,3);
n2 = coef(2,1)*x + coef(2,2)*y + coef(2,3);
n3 = coef(3,1)*x + coef(3,2)*y + coef(3,3);

% n = a * [x;y;1];
% n = n.double;

%First column is value, second column is df/dx, second column is df/dy
n = [n1.double; n2.double; n3.double];

% Jacobian
J = [n(1,2)*coord(1,1) n(1,3)*coord(2,1); n(2,2)*coord(1,2) n(2,3)*coord(2,2)];
disp(J);