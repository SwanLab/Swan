clear;

%coordinates of the quadrilater vertex (x1, y1, x2, y2, etc) (random)
coordElem = [ 1 1; 3 2; 5 4; 2 5];

N = zeros(8);

xiCoord = [ -1 1 1 -1 0 1 0 -1 ];
etaCoord = [ -1 -1 1 1 -1 0 1 0 ];

syms xi
syms eta

for i = 1 : 4

N(i) = N(i) + 0.25 * ( 1 + xiCoord(i) * xi ) * ( 1 + etaCoord(i) * eta ) * ( xiCoord(i) * xi + etaCoord(i) * eta - 1 );

end

for i = 5 : 7

N(i) = N(i) + 0.5 * ( 1 - xi^2 ) * ( 1 + etaCoord(i) * eta );

end

for i = 6 : 8

N(i) = N(i) + 0.5 * ( 1 + xiCoord(i) * xi ) * ( 1 - eta^2 );

end