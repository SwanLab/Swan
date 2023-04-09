clear;

%coordinates of the triangle vertex (x1, y1, x2, y2, etc) (random)
coordElem = [ 2 3; 6 2; 7 8; 1 7];

coef = zeros(2,4);

% matrix of coefficients
for i = 1 : 2 %coefficients of Nj *  parametric coordinates. Linear equation

coef(i,1) = 0.25 * ( coordElem(1,i) - coordElem(2,i) + coordElem(3,i) - coordElem(4,i) );
coef(i,2) = 0.25 * ( - coordElem(1,i) + coordElem(2,i) + coordElem(3,i) - coordElem(4,i) );
coef(i,3) = 0.25 * ( - coordElem(1,i) - coordElem(2,i) + coordElem(3,i) + coordElem(4,i) );
coef(i,4) = 0.25 * ( coordElem(1,i) + coordElem(2,i) + coordElem(3,i) + coordElem(4,i) );

end

%AD using my code
xi = ValDerForward(0,[1 0]);
eta = ValDerForward(0,[0 1]);

coord(1) = coef(1,1) * xi * eta + coef(1,2) * xi + coef(1,3) * eta + coef(1,4);
coord(2) = coef(2,1) * xi * eta + coef(2,2) * xi + coef(2,3) * eta + coef(2,4);

N = [ coord(1).double; coord(2).double ];

%Jacobian
J = transpose(N(:,2:end));

detJ = det(J);