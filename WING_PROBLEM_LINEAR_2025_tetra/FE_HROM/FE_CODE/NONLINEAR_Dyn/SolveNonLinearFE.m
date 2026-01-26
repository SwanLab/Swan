function[NODES_SNAP,GAUSS_SNAP,OPERfe,DATA,DATAOUT,NODESV_PROP,GAUSSV_PROP]  = ...
    SolveNonLinearFE(COOR,CN,TypeElement,TypeElementB, PROPMAT,  DOFr,dR,...
    Tnod,CNb,fNOD,Fpnt,typePROBLEM,DATA,CONNECTb,DOFm,Gb,MaterialType) ;

% Copy of SolveElastFE.m (adapted for nonlinear problems). JAHO,
% 17-Sept-2018
% -------------
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
%  .... Only valid for ELASTIC
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
% A) Global stiffness= matrix
% ------------------------------
COMP_TIME = tic ; 

DATA = DefaultField(DATA,'ComputeBstW',0) ; %    10-July-2019
 
disp('Computing B-matrices...')
 [wSTs,  XeALL,  Bst, BstW, wST,  DATA] = ...
    ComputeBmatrices(COOR,CN,TypeElement,DATA);

%------------------------

% B) External force vector due to body forces
% ------------------------------
disp('Computing   external force vector due to body forces (Fb)...')
[Fb, Nst,  DATA, wSTs_RHS, posgp_RHS] = BodyExternalForcesGlo(COOR,CN,TypeElement, fNOD,DATA) ;

% C)  External force vector due to   boundary tractions
% ------------------------------
disp('Computing  external force vector due to   boundary tractions ..')
Ftrac = TractionForcesGlo(COOR,CNb,TypeElementB,Fpnt,Tnod,CONNECTb,DATA);
ndim = size(COOR,2) ;

DATA = DefaultField(DATA,'PLOT_INPUT_FORCES',1) ; 

if DATA.PLOT_INPUT_FORCES == 1
    FinputNODES = Fb+Ftrac ; 
else
    FinputNODES = [] ; 
end
save(DATA.nameWORKSPACE,'FinputNODES','-append') ; 


%dbstop('109')
if DATA.CALCULATE_MASSMATRIX == 1
    
    % if  DATA.RECALCULATE_STIFFNESS == 1
    M= MassMatrix(DATA,Nst,wSTs_RHS,ndim) ;
    % else
    disp('Storing Mass Matrix...')
    save(DATA.nameWORKSPACE,'M','-append');
    disp('Done')
    %end
else
    M = [] ; 
end

% D) Solving for the vector of unknown displacements
disp('Solving...')
ASSEMBLY_INFO = [] ; 

[NODES_SNAP,GAUSS_SNAP,OPERfe,NODESV_PROP,GAUSSV_PROP,DATA] =...
    SolveNONLINEAR(Fb,Ftrac,dR,DOFr,COOR,CN,TypeElement,PROPMAT,typePROBLEM,...
    Bst,DOFm,Gb,DATA,wST,Nst,BstW,M,ndim,ASSEMBLY_INFO) ;
TOTAL_TIME_NONLINEAR_ANALYSIS = toc(COMP_TIME) ; 

disp([' TOTAL TIME NONLINEAR ANALYSIS = ',num2str(TOTAL_TIME_NONLINEAR_ANALYSIS)] )

save(DATA.nameWORKSPACE,'TOTAL_TIME_NONLINEAR_ANALYSIS','-append') ; 


% DATAOUT.strainGLO = strainGLO ;
% DATAOUT.stressGLO = stressGLO ;
%DATAOUT.wSTs = wSTs ;
%DATAOUT.MaterialType = MaterialType ;
%DATAOUT.d = d ;
%DATAOUT.Bst = Bst;
%stressGLO = DATAOUT.stress;
DATA_INPUT_FE = DATA ;
d = NODES_SNAP.U ; 
stressGLO = GAUSS_SNAP.stressST ; 
DOFr = OPERfe.DOFs ; 
posgp = DATA.posgp ; 
save(DATA.nameWORKSPACE,'d','stressGLO','DOFr','TypeElement','posgp','DATA_INPUT_FE','-append')
DATAOUT.nameWORKSPACE = DATA.nameWORKSPACE ;
