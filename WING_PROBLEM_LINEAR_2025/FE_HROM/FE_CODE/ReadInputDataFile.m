function [COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,...
    CONNECTb,DOFm,Gb,DATA,MaterialType,DOMAINVAR]  = ReadInputDataFile(FUNinput,DATA) ;


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
% 7. CONNECTb : Boundary elements connectivity
% 8. DOMAINVAR ---> Used only when meshes are created by tiling copies of
% single slice

% default inputs
% --------------
%dbstop('49')
if nargin == 1
    DATA = [] ; 
end

[DATA DOFm Gb MaterialType]= SetDefaultInputs(DATA) ; 

disp('Reading input data...')

time1 = tic ;

NAME_INPUT_DATA = FUNinput.NAME ;
INPUTS_LOC =  FUNinput.INPUTS ;
FUNinput = DefaultField(FUNinput,'NAME_DATA_INPUT',''); % JAHO, 4-Apr-2019
DATA = DefaultField(DATA,'INPUTDATAfile',FUNinput.NAME_DATA_INPUT);  % JAHO, 21-Mar-2019
DATA = DefaultField(DATA,'REACTIONS_RESULTANTS_CALCULATE',[]);  % JAHO, 8-July-2019

 
 
[COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...  
    Tnod,CNb,fNOD,Fpnt,NameFileMesh,typePROBLEM,celasgloINV,CONNECTb,DOFm,Gb,...
    DATA,MaterialType] = feval(NAME_INPUT_DATA,INPUTS_LOC,DATA) ;

DATA = DefaultField(DATA,'DOMAINVAR',[]) ; 
DOMAINVAR= DATA.DOMAINVAR ; 
DATA.DOMAINVAR = [] ; 

% DATA = DefaultField(DATA,'DOMAINVAR',[]) ;
%  
% DOMAINVAR = DATA.DOMAINVAR  ; 
% DATA.DOMAINVAR = [] ; 


% %% RENUMERING CONNECTIVITY MATRIX (to ensure small bandwidth of both BB and BBnw)
% [~,IndicesRenumberingElements]  = sort(CN(:,1)) ;
% CN = CN(IndicesRenumberingElements,:) ;
% celasglo = celasglo(:,:,IndicesRenumberingElements) ;
% if ~isempty(MaterialType)
%     MaterialType = MaterialType(IndicesRenumberingElements) ;
% end

 
 
disp(['Done (in ',num2str(toc(time1)),' s)']);

