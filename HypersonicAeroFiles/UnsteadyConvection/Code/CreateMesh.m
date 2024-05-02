function [X,T]=CreateMesh(x1,x2,y1,y2,nx,ny)
% [X,T] = CreateMesh(x1,x2,y1,y2,nx,ny)
% Topology for a structured mesh of bilinear quadrilaterals
% in the rectangular domain [x1,x2]x[y1,y2]
% nx,ny: number of elements in eavh direction
%


X = zeros((nx+1)*(ny+1),2);
T = zeros(nx*ny,4);

hx = (x2-x1)/nx;
hy = (y2-y1)/ny;

xs = [x1:hx:x2]';
vect_ones = ones(nx+1,1);

% nodes coordinates
for i=1:ny+1
    ys = ((i-1)*hy+y1)*vect_ones;
    posi = [(i-1)*(nx+1)+1:i*(nx+1)]; 
    X(posi,:)=[xs,ys];
end

% Connectivities
for i=1:ny
    for j=1:nx
        ielem = (i-1)*nx+j;
        inode = (i-1)*(nx+1)+j;
        T(ielem,:) = [ inode inode+1 inode+(nx+2) inode+(nx+1)];
    end
end
