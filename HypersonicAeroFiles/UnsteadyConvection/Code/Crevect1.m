function f = Crevect1 (X,T,pospg,wpg,N,Nxi,Neta);
% f = CreVect1 (X,T,pospg,wpg,N,Nxi,Neta);
% Computation of vector f ohbtained by discretizing (w,s) 
% where s is definie in SourceTerm
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
    % Element Matrix
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
        % Contribution at the element vector
        aux = Isopar(Xe,N_igaus);
        f_igaus = SourceTerm(aux); 
        fe = fe + N_igaus'*f_igaus*dvolu; 
    end   
    % Assembly of the element vector
    f(Te) = f(Te) + fe; 
    clear fe; 
end 

