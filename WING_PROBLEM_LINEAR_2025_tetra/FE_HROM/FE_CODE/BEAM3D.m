% Inputs example assigment 2 
% ----------------------------
%%%%%%%%%%%%%%%%%
 % ---------------
% 1.  Finite element mesh:  COORDINATES AND CONNECTIVITIES for both the volume domain and the boundary domain
% OUTPUT: COOR,CN,TypeElement,CONNECTb,TypeElementB
%NameFileMesh = 'beam2elem.msh'; % Name of the file containing the mesh information (Generated with GID)
%NameFileMesh = 'mesh3D10e3.msh' ;
NameFileMesh = 'mesh3D10e5.msh';

[COOR,CN,TypeElement,CONNECTb,TypeElementB]=...
    ReadMeshFile(NameFileMesh)  ;

nnode = size(COOR,1) ;% Number of nodes 
ndim = size(COOR,2); % Number of spatial dimensions (ndim=2 for 2D problems)
nelem = size(CN,1) ; % Number of elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. MATERIAL PROPERTIES: output celasglo   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
typePROBLEM = 'pstress';  %'pstress'/'pstrain'/'3D';  Plane stress/ plane strain problem 
if ndim==2 
    nstrain = 3; 
else
    nstrain = 6 ; 
    typePROBLEM ='3D' ;
end
celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
celasgloINV = zeros(6,6,nelem) ;  % Global array of compliance matrices (3D)

% Elastic properties (isotropic)
E = 3e7  ; %  Pa, Young's modulus
nu = 0.3; % Poisson's coefficient
% Compliance matrix for an isotropic materials (with all entries, 3D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = E/2/(1+nu) ; 
celasINV3D = [1/E  -nu/E -nu/E  0 0 0  
           -nu/E  1/E  -nu/E  0 0 0
           -nu/E  -nu/E  1/E  0 0 0           
            0      0    0  1/G   0 0
            0      0      0  0 1/G 0  
            0      0      0  0  0  1/G] ; 
switch typePROBLEM
    case 'pstrain'
        celas3D = inv(celasINV3D)  ; 
        rowcol = [1 2 6] ; 
        celas = celas3D(rowcol,rowcol) ; 
    case 'pstress' 
        rowcol = [1 2 6] ; 
        celasINV = celasINV3D(rowcol,rowcol) ; 
        celas = inv(celasINV) ; 
    case '3D'
        celas = inv(celasINV3D)  ; 
end 
for e=1:nelem 
    celasglo(:,:,e) = celas ; 
    celasgloINV(:,:,e) = celasINV3D ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Dirichlet (essential) boundary conditions, OUTPUT: dR and rdof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rnod = {} ; uPRES={}  ;
% 1) List of nodes at which displacement is prescribed  in DIRECTION  i = 1
% All nodes pertaining to plane z = 0 
% Boundary nodes 
BoundaryNodes= unique(CONNECTb(:)) ;
xmin = min(COOR(BoundaryNodes,1)) ; xmin = xmin(1) ;
idim=1 ; 
rnodBASEloc = find(abs(COOR(BoundaryNodes,1)-xmin)<1e-10) ; 
rnodBASE = BoundaryNodes(rnodBASEloc) ; 
rnod{idim} =rnodBASE;  
% Vector of prescribed displacements
displ1 = 0 ; 
uPRES{idim} = displ1*ones(size(rnod{idim})) ; 
% 2) List of nodes at which displacement is prescribed  in DIRECTION  i = 2
idim = 2;  
rnod{idim} =rnodBASE ;  
% Vector of prescribed displacements
displ2 = 0 ; 
uPRES{idim} = displ2*ones(size(rnod{idim})) ; 
% 2) List of nodes at which displacement is prescribed  in DIRECTION  i = 3
idim = 3;  
rnod{idim} =rnodBASE ;  
% Vector of prescribed displacements
displ3 = 0 ; 
uPRES{idim} = displ3*ones(size(rnod{idim})) ; 
%%%% Set of restricted degrees of freedom and vector of prescribed
%%%% displacements (dR)
DOFr = [] ; dR = [] ; 
for idim = 1:ndim 
    DOFr = [DOFr ; (rnod{idim}-1)*ndim+idim]; 
    dR = [dR ; uPRES{idim}]; 
end

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Neumann (natural) boundary conditions : OUTPUT: Tnod, CNb, Fnod  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POINT LOADS 
% -----------
Fpnt =zeros(ndim*nnode,1) ; % There are no point loads
% DISTRIBUTED LOADS
% ------------------------
CNb ={} ; Tnod={} ;
% 1) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 1
%                      %%%%%%%%
idim = 1; 
CNb{idim} = [] ;   % In this case, there are no distributed loads in the idim=1 direction
Tnod{idim} = [] ;  
% ----
% 2) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 2
%                      %%%%%%%%
idim = 2; 
% Top surface 
ymax = max(COOR(BoundaryNodes,2)) ; ymax = ymax(1) ; 
nodTOPloc = find(abs(COOR(BoundaryNodes,2)-ymax)<1e-10) ; 
NODESb = BoundaryNodes(nodTOPloc) ; 
% Connectivity matrix for boundary elements (Neumann boundary)
CNb{idim} = ElemBnd(CONNECTb,NODESb) ; 
% Value of the distributed load at the nodes of each boundary element 
ty = -20 ;  % N/m
Tnod{idim} = ty*ones(size(CNb{idim})) ; 
% 3) List of boundary nodes whose corresponding boundary elements have distributed loads
% along  direction --> idim = 3
%       
idim = 3; 
CNb{idim} = [] ;   % In this case, there are no distributed loads in the idim=1 direction
Tnod{idim} = [] ;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Heat source
%% 
f = 0 ; % kN/m3
fNOD =  zeros(nnode*ndim,1) ; 
fNOD(2:3:3*nnode) = f ; 
