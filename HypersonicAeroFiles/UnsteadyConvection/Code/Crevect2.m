function f = Crevect2 (X,T,Conv,pospg,wpg,N,Nxi,Neta);
% f = Crevect2 (X,T,Conv,pospg,wpg,N,Nxi,Neta);
% Computation of vector f ohbtained by discretizing (a·grad(w),s) 
% where s is defined in SourceTerm
% 

% Total number of elements and number of nodes in each element
[numel,nen] = size(T); 
% Total number of nodes
numnp = size(X,1); 
% Number of Gauss points in one element
ngaus = size(pospg,1); 

% Allocation 
f = zeros(numnp,1); 
 
% Loop on elements
for ielem = 1:numel 
    % Te: global number of nodes on the current element
    Te = T(ielem,:); 
    % Xe: coordinates of the nodes in Te
    Xe = X(Te,:); 
    % Conve: velocity field on Te
    Conve = Conv(Te,:);    
    % Element vector
    fe = zeros(nen,1); 
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
        % Derivatives of shape functions on global coordinates
        res = Jacob\[Nxi_igaus;Neta_igaus]; 
        Nx = res(1,:); 
        Ny = res(2,:); 
        % Contribution at the element vector
        a = Isopar(Conve,N_igaus); 
        aux = Isopar(Xe,N_igaus);
        f_igaus = SourceTerm(aux); 
        fe = fe + (a(1)*Nx + a(2)*Ny)'*f_igaus*dvolu; 
    end   
    % Assembly of the element vector
    f(Te) = f(Te) + fe; 
    clear fe; 
end 

