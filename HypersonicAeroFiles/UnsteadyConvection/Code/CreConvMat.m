function C = CreConvMat (X,T,Conv,pospg,wpg,N,Nxi,Neta);
% C = CreConvMat (X,T,Conv,pospg,wpg,N,Nxi,Neta);
% Computation of the convection matrix C 
% obtained by discretizing (a·grad(w), u)
%


% Total number of elements and number of nodes in each element
[numel,nen] = size(T); 
% Total number of nodes
numnp = size(X,1); 
% Number of Gauss points in one element
ngaus = size(pospg,1); 

% Allocation 
C = zeros(numnp,numnp); 
 
% Loop on elements
for ielem = 1:numel 
    % Te: global number of nodes on the current element
    Te = T(ielem,:); 
    % Xe: coordinates of the nodes in Te
    Xe = X(Te,:); 
    % Conve: velocity field on Te
    Conve = Conv(Te,:);
    % Element Matrix
    Ce = zeros(nen,nen); 
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
        Ce = Ce + (a(1)*Nx+a(2)*Ny)'*N_igaus*dvolu; 
    end   
    % Assembly of the element matrix
    C(Te,Te) = C(Te,Te) + Ce; 
    clear Ce; 
end 

C = sparse(C); 