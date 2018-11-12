% PERIODIC BOUNDARY COND
% creation of a list containing the couples that define the periodicity
% 1) list of nodes on each side (two vertical, two horizontal)
% 2) sort each list based on the corresponding coordinate
% 3) generation of the couples and store them in 'element.pnods'
% 
% remark: 
% [Y,I] = SORT(X,DIM,MODE) also returns an index matrix I.
% If X is a vector, then Y = X(I).  

% nodes in the left-vertical side, without the corners
href = 0.025;
h=href; L=[];Y=[];
for i=1:size(coordinates,1)
    if (coordinates(i,1)<h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5 )
        L = [L i];
        Y = [Y coordinates(i,2)];
    end
end
[Y1,I] = sort(Y);
V1 = L(I);

% nodes in the right-vertical side, without the corners
h=href; L=[];Y=[];
for i=1:size(coordinates,1)
    if (coordinates(i,1)>1-h/5 && coordinates(i,2)>h/5 && coordinates(i,2)<1-h/5)
        L = [L i];
        Y = [Y coordinates(i,2)];
    end
end
[Y1,I] = sort(Y);
V2 = L(I);

% nodes in the bottom-horizontal side, without the corners
h=href; L=[];X=[];
for i=1:size(coordinates,1)
    if (coordinates(i,2)<h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5 )
        L = [L i];
        X = [X coordinates(i,1)];
    end
end
[X1,I] = sort(X);
H1 = L(I);

% nodes in the top-horizontal side, without the corners
h=href; L=[];X=[];
for i=1:size(coordinates,1)
    if (coordinates(i,2)>1-h/5 && coordinates(i,1)>h/5 && coordinates(i,1)<1-h/5)
        L = [L i];
        X = [X coordinates(i,1)];
    end
end
[X1,I] = sort(X);
H2 = L(I);
element.pnods = [V1 H1; V2 H2]; % lista de nodos 

element.penalty = 0;
% element.lglib = [2*1-1    2*3277-1  2*1     2*3277 ;  ...
%                  2*4225-1 2*3278-1  2*4225  2*3278];
% element.l2penalty = [   1  3277    1 3278;
%                      3278  4225 3277 4225];
% element.penalty = 10000000*1*0.015625^2;

%  1    3277
%
% 3278  4225

