clear;

% symbolic variables
syms x
syms y

%coordinates of the triangle vertex (x1,y1; x2,y2; etc) (random)
coord = [2 2; 7 1; 5 4];

% xi and eta equation, found on article Chapter 4 Finite Element
% Approximation page 6
xi = (det([1 x y; 1 coord(1,1) coord(1,2); 1 coord(3,1) coord(3,2)]))/(det([1 coord(2,1) coord(2,2); 1 coord(1,1) coord(1,2); 1 coord(3,1) coord(3,2)]));
eta = (det([1 x y; 1 coord(1,1) coord(1,2); 1 coord(2,1) coord(2,2)]))/(det([1 coord(3,1) coord(3,2); 1 coord(1,1) coord(1,2); 1 coord(2,1) coord(2,2)]));

N = [xi eta 1-xi-eta];

%transform symbolic equation into "normal" equation
a = [double(coeffs(N(1))); double(coeffs(N(2))); double(coeffs(N(3)))];

% AD using my code.
x = ValDerForward(1,[1 0]);
y = ValDerForward(1,[0 1]);

n1 = a(1,1)*x + a(1,2)*y + a(1,3);
n2 = a(2,1)*x + a(2,2)*y + a(2,3);
n3 = a(3,1)*x + a(3,2)*y + a(3,3);

% n = a * [x;y;1];
% n = n.double;

%First column is value, second column is df/dx, second column is df/dy
n = [n1.double; n2.double; n3.double];

% Jacobian
J = [n(1,2)*coord(1,1) n(1,3)*coord(2,1); n(2,2)*coord(1,2) n(2,3)*coord(2,2)];
%
% disp(J);