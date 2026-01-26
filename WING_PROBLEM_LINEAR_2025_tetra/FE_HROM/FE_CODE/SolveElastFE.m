function[d, strainGLOgid, stressGLOgid,   React, posgp, CN, MaterialType, DATAOUT,Fnodes]  = ...
    SolveElastFE(COOR,CN,TypeElement,TypeElementB, celasglo,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,typePROBLEM,celasgloINV,DATA,CONNECTb,DOFm,Gb,MaterialType) ;

%%% This function returns the (nnode*ndim x 1) vector of nodal displacements (d),
%%% as well as the arrays containing  the stresses (stressGLO) and strains (strainGLO) at all gauss
%%% points
% % INPUTS
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
%  fNOD: Vector containing the nodal values of the heat source function (nnode*ndime x1 )%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Miscellaneous input DATA  -->
% DATA.VECTcode = 1 --> VEctorized code
%
d=[]; strainGLO=[] ; stressGLO=[] ;posgp=[] ;
% A) Global stiffness= matrix
% ------------------------------
disp('Computing stiffness matrix K ...')
[K, CN, wSTs,  XeALL, Cglo, Bst, wST, MaterialType, DATA] = ...
    StiffnessMatrixGlobal(COOR,CN,TypeElement,celasglo,DATA,MaterialType);


%------------------------

% B) External force vector due to body forces
% ------------------------------
aaaa = tic;
disp('************************************************************+')
disp('Computing   external force vector due to body forces (Fb)...')
disp('************************************************************+')

[Fb, Nst,  DATA, wSTs_RHS, posgp_RHS] = BodyExternalForcesGlo(COOR,CN,TypeElement, fNOD,DATA) ;
aaaa =toc(aaaa) ;
disp([' DONE (in ',num2str(aaaa),' seconds)']) ;
disp('************************************************************+')

% C)  External force vector due to   boundary tractions
% ------------------------------
disp('************************************************************+')

aaaa = tic;
disp('Computing  external force vector due to   boundary tractions ..')
disp('************************************************************+')
aaaa =toc(aaaa) ;

Ftrac = TractionForcesGlo(COOR,CNb,TypeElementB,Fpnt,Tnod,CONNECTb,DATA);

disp([' DONE (in ',num2str(aaaa),' seconds)']) ;
disp('************************************************************+')

Fnodes = Fb + Ftrac ;


%dbstop('109')
if DATA.CALCULATE_MASSMATRIX == 1
    ndim = size(COOR,2) ;
    % if  DATA.RECALCULATE_STIFFNESS == 1
    M= MassMatrix(DATA,Nst,wSTs_RHS,ndim) ;
    densGLO = DATA.densGLO;
    % else
    disp('Storing Mass Matrix...')
    save(DATA.nameWORKSPACE,'M','densGLO','-append');
    disp('Done')
    %end
end

% D) Solving for the vector of unknown displacements
DATA = DefaultField(DATA,'MINIMIZE_MEMORY_REQUIREMENTS_IN_SOLVING',3) ;

disp('Solving...')
if DATA.NOCALCULATE_DISPLACEMENTS == 0
    [d strainGLOgid stressGLOgid  React posgp DATAOUT] =...
        SolveELAS(K,Fb,Ftrac,dR,DOFr,COOR,CN,TypeElement,celasglo,typePROBLEM,celasgloINV,...
        Cglo,Bst,DOFm,Gb,DATA,wST,Nst) ;
else
    DATAOUT.stress = [] ;     strainGLOgid = [] ;     stressGLOgid = [] ;
    d = [] ; React = [] ;
end

% DATAOUT.strainGLO = strainGLO ;
% DATAOUT.stressGLO = stressGLO ;
DATAOUT.wSTs = wSTs ;
DATAOUT.MaterialType = MaterialType ;
DATAOUT.d = d ;

DATA = DefaultField(DATA,'dispMACRO',[] ) ; % MAcro-displacement (homogenization problems)

if ~isempty(DATA.dispMACRO)
    DATAOUT.dispMACRO=  DATA.dispMACRO ;
end

%DATAOUT.Bst = Bst;
DATAOUT = DefaultField(DATAOUT,'stress',[]) ;
stressGLO = DATAOUT.stress;
DATA_INPUT_FE = DATA ;
if DATA.STORE_STIFFNESS ~=0
    save(DATA.nameWORKSPACE,'d','stressGLO','DOFr','TypeElement','posgp',...
        'DATA_INPUT_FE','Fnodes','React','DOFm','Gb','-append')
end

DATAOUT.nameWORKSPACE = DATA.nameWORKSPACE ;
