function M = CreMassMat (X,T,pospg,wpg,N,Nxi,Neta);
% M = CreMassMat (X,T,pospg,wpg,N,Nxi,Neta);
% Computation of the mass matrix M
% ohbtained by discretizing (w, u)
%


% Total number of elements and number of nodes in each element
[numel,nen] = size(T); 
% Total number of nodes
numnp = size(X,1); 
% Number of Gauss points in one element
ngaus = size(pospg,1); 

% Allocation 
M = zeros(numnp,numnp); 
 
% Loop on the elements
for ielem = 1:numel 
    % Te: global number of nodes on the current element
    Te = T(ielem,:); 
    % Xe: coordinates of the nodes in Te
    Xe = X(Te,:); 
    % Element Matrix
    Me = zeros(nen,nen); 
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
        % Contribution at the element matrix
      Me = Me + N_igaus'*N_igaus*dvolu; 
    end   
    % Assembly of the element matrix
    M(Te,Te) = M(Te,Te) + Me; 
    clear Me; 
end 

M = sparse(M); 

