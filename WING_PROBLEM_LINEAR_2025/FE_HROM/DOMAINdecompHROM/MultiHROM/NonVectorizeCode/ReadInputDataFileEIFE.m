function [COOR,CN,TypeElement,TypeElementB,  DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,DATA,PROPMAT,MaterialType]  = ReadInputDataFileEIFE(NAME_INPUT_DATA,DATA) ; 
% EIFE METHOD, READING INPUTS, 9-mARCH-2023, jaho

% OUTPUTS 
% --------------
% 1. Finite element mesh 
% -------------------
% COOR: Coordinate matrix (nnode x ndim)
% CN: Connectivity matrix (nelem x nnodeE)
% TypeElement: Type of finite element (quadrilateral,...)
% TypeElementB: Type of boundary finite element (linear...)
% -----------
% 2. Material
% -----------
%  celasglo (nstrain x nstrain x nelem)  % Array of elasticity matrices
%  celasgloINV (6 x 6 x nelem)  % Array of compliance matrices (3D)
% -------------------------
% 3. Dirichlet (Essential) Boundary Condition s
% --------------------------------------------
%  DOFr --> Set of Global Degrees of Freedom with prescribed displacements 
%  dR   --> Vector of prescribed displacements  (size(DOFr) = size(dR))
% ---------------------------------------
% 4. Neumann (natural) Boundary Conditions
% -------------------------------------------
% DISTRIBUTED LOADS
% -------------------
%  CNb: Cell array in which the cell CNb{idim} contains the connectivity matrix for the boundary elements
%  of the traction boundaries in the idim direction
%  Tnod: Cell array in which the entry Tnod{idim} features the vector with the prescribed traction at
%   the nodes specified in CNb{idim}    (note that size(CNb{idim}) = size(Tnod{idim}))
%  each cell of 
%  POINT LOADS
% --------------------------------
%  Fpnt  (nnode*ndime x 1):  Vector containing point forces applied on the
%  nodes of the discretization
% ----------------------------------
% 5. Body force
% ---------------
%  fNOD: Vector containing the nodal values of the heat source function (nnode*ndime x1 )
%% 6. Type problem  (plain stress/strain/3D)
% typePROBLEM = 'pstress'/'pstrain'; 
 disp('Reading input data...')
 celasgloINV = [] ; densglo = [] ; 
 %DATA = [] ; 
eval(NAME_INPUT_DATA) ; 
