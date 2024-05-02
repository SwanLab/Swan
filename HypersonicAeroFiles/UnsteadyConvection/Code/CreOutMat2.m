function C = CreOutMat2 (X,T,Conv,elem_out,side);
% C = CreOutMat2 (X,T,Conv,elem_out,side);
% Computation of the matrix M obtained by discretizing 
% the term a_n·(w,a·grad(u)) over the boundary
%


% Total number of elements 
numel = length(elem_out); 
% Total number of nodes
numnp = size(X,1); 

% number of Gauss points (numerical quadrature on the boundary)
ngaus = 2; 
% number of shape functions on each element 
nfunc = 4;
% Gauss points and weigths for the reference element [-1,1]
[pospg,wpg] = Quadrature_cont (ngaus,side);
% Shape functions for the linear 1D element
[N,Nxi,Neta] = ShapeFunc(pospg);

% Computation of unitary tangent vector taking into account that
% boundaries are parallel to one of the coordinate axes
Te = T(elem_out(1),:);  
Xe = X(Te(side),:);   
x1 = Xe(1,:); x2 = Xe(ngaus,:);
vtan = x2-x1;
h = norm(x2-x1);
h_2 = h/2;
vnorm = [vtan(2), - vtan(1)]/h;  

% Allocation 
C = zeros(numnp,numnp); 
 
% Loop on the elements
for ielem = 1:numel 
    % Te: global number of nodes on the current element
    Te = T(elem_out(ielem),:); 
    % Xe: coordinates of the nodes in Te
    Xe = X(Te,:); 
    % Conve: velocity field on Te
    Conve = Conv(Te,:);
    % Element Matrix
    Ce = zeros(nfunc,nfunc); 
    % Loop on Gauss points (numerical quadrature)
    for igaus = 1:ngaus
        % Shape functions on Gauss point igaus
        N_igaus    = N(igaus,:);
        Nxi_igaus  = Nxi(igaus,:);
		Neta_igaus = Neta(igaus,:);
        %
        dvolu = wpg(igaus)*h_2;
        % Derivatives of the shape functions on global coordinates
        Nx = Nxi_igaus/h_2;
        Ny = Neta_igaus/h_2;
        % velocity on the Gauss point
        a = Isopar(Conve,N_igaus); 
        % normal component of the velocity o the Gauss point
        an = a*vnorm';
        % Contribution to element matrix
        Ce = Ce + an*N_igaus'*(a(1)*Nx+a(2)*Ny)*dvolu;
    end   
    % Assembly of the element matrix
    C(Te,Te) = C(Te,Te) + Ce; 
    clear Ce; 
end 

C = sparse(C); 

