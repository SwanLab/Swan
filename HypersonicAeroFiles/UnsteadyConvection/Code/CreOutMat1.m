function M = CreOutMat1 (X,T,Conv,elem_out,side);
% M = CreOutMat1 (X,T,Conv,elem_out,side);
% Computation of the matrix M obtained by discretizing 
% the term a_n·(w,u) over the outflow boundary
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
M = zeros(numnp,numnp); 
 
% Loop on elements
for ielem = 1:numel 
    % Te: global number of nodes on the current element
    Te = T(elem_out(ielem),:); 
    % Xe: coordinates of the nodes in Te
    Xe = X(Te,:); 
    % Conve: velocity field on Te
    Conve = Conv(Te,:);
    % Element Matrix
    Me = zeros(nfunc,nfunc); 
    % Loop on Gauss points (numerical quadrature)
    for igaus = 1:ngaus
        % Shape function on Gauss point igaus
        N_igaus = N(igaus,:);
        %
        dvolu = wpg(igaus)*h_2;
        % velocity on the Gauss point
        a = Isopar(Conve,N_igaus); 
        % normal component of the velocity on the Gauss point
        an = a*vnorm';
        % Contribution to element matrix
        Me = Me + an*N_igaus'*N_igaus*dvolu;
    end   
    % Assembly of the element matrix
    M(Te,Te) = M(Te,Te) + Me; 
    clear Me; 
end

M = sparse(M); 

