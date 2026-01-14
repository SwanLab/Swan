% Inputs example assigment 2 
% ----------------------------


%%%%%%%%%%%%%%%%%
 % ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB
NameFileMesh = 'malla1.msh'; % Name of the file containing the mesh information (Generated with GID)
[COOR,CN,TypeElement,CONNECTb,TypeElementB]=...
    ReadMeshFile(NameFileMesh)  ;
nnode = size(COOR,1) ;% Number of nodes 

% 2. MATERIAL PROPERTIES: output ConductMglo  
%-----------------------
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements
ConductMglo = zeros(ndim,ndim,nelem) ; 
% Conductivity matrix (isotropic)
kappa =  5  ; %  W/ºC
ConductM = kappa*eye(ndim) ; % eye = IDENTITY MATRIx
for e=1:nelem 
    ConductMglo(:,:,e) = ConductM ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rnod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of nodes at which temperature is prescribed
%   Copy here the list of nodes from GID
rnod =[1         0.0 0.0 0.0
    2         0.0 0.25 0.0
    3         0.0 0.5 0.0
    4         0.5 0.125 0.0
    6         0.0 0.75 0.0
    9         0.0 1.0 0.0
    10         1.0 0.25 0.0
    16         1.5 0.375 0.0
    21         2.0 0.5 0.0];  %

rnod =rnod(:,1) ;  
% Prescribed temperature 
temp_DAB = 0 ; % Degrees Celsius 
% Vector of prescribed temperatures
dR = temp_DAB*ones(size(rnod)) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: qFLUXglo, CNb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of boundary nodes whose boundary elements have prescribed heat flux
% (different from zero !)
% Edge DC 
% Copy here the list of nodes from GID
NODESb = [9         0.0 1.0 0.0
12         0.5 1.0 0.0
15         1.0 1.0 0.0
20         1.5 1.0 0.0
25         2.0 1.0 0.0] ; 
%
NODESb = NODESb(:,1) ; 
% prescribed flux (constant)
qBAR = 20 ; % W/m
% Initialization  of the global vector of nodal heat flux 
qFLUXglo = zeros(nnode,1) ; 
% Assign prescribed values
qFLUXglo(NODESb) = qBAR ; 
%%%%
% Connectivity matrix for boundary elements (Neumann boundary)
CNb = ElemBnd(CONNECTb,NODESb) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Heat source
%%
f = 6 ; %   W/m2
% Global vector of heat sources (constant)
fNOD = f*ones(nnode,1) ; 
