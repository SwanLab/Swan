function CNb =  IdentifyBoundaryMesh(CN)
% Given Connectivity matrix CN of a 2D mesh, this function returns the connectivities of the boundary nodes
%(triangle and quadrilateral meshes)
% Joaquin A. Hernandez Ortega (JAHO), 8-April-2018
if nargin == 0
    load('tmp.mat')
end
nnodeE = size(CN,2);  % Number of nodes 

if nnodeE == 4 
    % QUADRILATERAL 
    % ----------------
    % Find all edges in the mesh
    E = [CN(:,1), CN(:,2)     % Edge 1
         CN(:,2), CN(:,3)    % Edge 2
         CN(:,3), CN(:,4)     % Edge 3
         CN(:,4), CN(:,1) ] ; % Edge 4
     
elseif nnodeE == 3 
      E = [CN(:,1), CN(:,2)     % Edge 1
         CN(:,2), CN(:,3)    % Edge 2
         CN(:,3), CN(:,1)   % Edge 3
          ] ; % Edge 4
end
% Find all edges in mesh, note internal edges are repeated
E = sort(E')';  % 
% determine uniqueness of edges
[u,m,n] = unique(E,'rows');
% determine counts for each unique edge
counts = accumarray(n(:), 1);
% extract edges that only occurred once
CNb = u(counts==1,:);

