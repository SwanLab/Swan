function K = CreStiffMat (X,T,Conv,pospg,wpg,N,Nxi,Neta);
% K = CreStiffMat (X,T,Conv,pospg,wpg,N,Nxi,Neta);
% Computation of the stiffness matrix K
% obtained by discretizing (a·grad(w), a·grad(u))
%


% Total number of elements and number of nodes in each element
[numel,nen] = size(T); 
% Total number of nodes
numnp = size(X,1); 
% Number of Gauss points in one element
ngaus = size(pospg,1); 

% Allocation 
K = zeros(numnp,numnp); 
 
% Loop on elements
for ielem = 1:numel 
    % Te: global number of nodes on the current element
    Te = T(ielem,:); 
    % Xe: coordinates of the nodes in Te
    Xe = X(Te,:); 
    % Conve: velocity field on Te
    Conve = Conv(Te,:);
    % Element Matrix
    Ke = zeros(nen,nen); 
    % Loop on Gauss points (numerical quadrature)
    for igaus = 1:ngaus
        % Shape functions on Gauss point igaus
        N_igaus    = N(igaus,:);     
        Nxi_igaus  = Nxi(igaus,:);   
        Neta_igaus = Neta(igaus,:);
        % Jacobian matrix on the Gauss point
        Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))   
                Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
        %
        dvolu = wpg(igaus)*det(Jacob); 
        % Derivatives of shape functions in global coordinates
        res = Jacob\[Nxi_igaus;Neta_igaus]; 
        Nx = res(1,:); 
        Ny = res(2,:); 
        % Contribution at the element matrix
        a = Isopar(Conve,N_igaus); 
        Ke = Ke + (a(1)*Nx+a(2)*Ny)'*(a(1)*Nx+a(2)*Ny)*dvolu; 
    end   
    % Assembly of the element matrix
    K(Te,Te) = K(Te,Te) + Ke; 
    clear Ke; 
end 

K = sparse(K); 