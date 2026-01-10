function [D_BmatIst,Nmat,WEIGHTSinteg,TRANSF_COORD,CN,DATA,D_KcLINloc,CN_SUPPORT,QrotINI,lambdaLEN,D_Wecm,INDEXsparseROTmat,...
    LboolCall,XcALL,INDEXsparseFINT,D_PdownsRBlROT,XcALLloc,D_AspMAT,INDEXsparseROTelem_byROWS] = B_N_matricesEIFEbubCOROT_LRss(COOR,CN, PROPMAT,MaterialType,DATA,TypeElement,...
    CN_SUPPORT)
%--------------------------------------------------------------------------
%  B_N_matricesEIFEbubCOROT_LRss
%
%  Constructs interscale matrices (B-matrix, N-matrix), weights, and rotation
%  operators for each EIFEM element under the small-strains/large-rotations 
%  regime. This function implements a **corotational extension** of the 
%  original EIFEM approach, allowing domain-wise consistency of rotational 
%  transformations as described in Sections 4, 6, and 12.8 of the document
%  "EIFEM_largeROTfinal.pdf".
%
%  INPUTS:
%    - COOR, CN: nodal coordinates and connectivity matrix
%    - PROPMAT: structure array with trained EIFE operators and integration rules
%    - MaterialType: vector assigning material models to elements
%    - TypeElement: element type (e.g., Q4, T3, etc.)
%    - CN_SUPPORT (optional): for uncoupled meshes (multiphysics or post-process)
%
%  OUTPUTS:
%    - D_BmatIst: Block-diagonal interscale strain operators
%    - Nmat: Global matrix interpolating coarse-scale displacements
%    - WEIGHTSinteg: Integration weights for internal/body forces
%    - TRANSF_COORD: Element-wise coordinate transformations, including initial Q
%    - D_KcLINloc: Local linear stiffness matrices from training
%    - QrotINI: Initial rotation matrices (per element)
%    - lambdaLEN: Scaling lengths for each element (used in XcALL)
%    - D_Wecm: Diagonal matrix of weights for internal forces
%    - D_PdownsRBlROT: Downscaling operator for rotation-induced displacements
%    - D_AspMAT: Spin operator matrices for large-rotation consistency
%    - INDEXsparse*: Index structures for efficient assembly of:
%        - Internal forces (FINT)
%        - Extended rotation matrices (global and per row in 3D)
%    - LboolCall: Boolean operator for global-to-element mapping
%    - XcALL, XcALLloc: Global and local reference axes
%    - DATA: enriched with mesh and ECM indexing
%
%  The routine:
%   - Chooses the best parent domain via `ParentDomainSearchEIFEcorot`
%   - Builds rigid-body-consistent B and N matrices with associated weights
%   - Prepares all quantities needed for consistent assembly in nonlinear regimes
%   - Computes support structures for corotational updates (e.g. Spin, Downscaling)
%   - Handles variable number of bubble DOFs across different element types
%
%  For each element, the output matrices reflect local transformations into 
%  reference domains (Sec. 6.3), and are rescaled via `lambdaLEN` and rotated 
%  using `QrotINI`. The global reference axes (`XcALL`) and spin matrices 
%  are used in the variational formulation of Section 12.
%
%  The output of this function is essential for:
%    - Consistent computation of internal forces and tangent matrices
%    - Tracking of element-level rotations and their variational coupling
%    - Efficient block-diagonal assembly in both training and online phases
%
%  Author:
%    Joaquín A. Hernández Ortega, UPC/CIMNE
%    Balmes 185, Barcelona
%    Version: 05-Feb-2025
%    Comments by ChatGPT4, 13-May-2025
%
%  Related scripts and implementations:
%    - 109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
%    - Bmat_weights_EIFEbubCOROT.m, Nmat_weights_EIFEbubCOROT.m
%
%--------------------------------------------------------------------------




% EIFE METHOD, DETERMINATION OF B and N matrices  --------------
% Adaptation of B_N_matricesEIFEbub.m
% CO-ROTATIONAL APPROACH,
% FORMULATION small strains/Large rotations (domain-wise)
% /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% 5-Feb-2O24, BALMES 185, BARCELONA
% --------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
nelem = size(CN,1); % Number of coarse-scale elements
TRANSF_COORD = cell(nelem,1) ; % INFO ABOUT COORDINATE TRANSFORMATION; to be used in the reconstruction phase
D_BmatIst = cell(nelem,1) ; % Interscale B-matrix Matrix that maps coarse-scale DOFS onto fine-scale strains at the ECM points of each element
Nmat  = cell(nelem,1) ; % Matrix that maps coarse-scale DOFS onto fine-scale displacements at the ECM points of each element
D_KcLINloc = cell(nelem,1) ; %  Linear stiffness matrix %/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
QrotINI  = cell(nelem,1) ;  % Initial rotation matrices of each EIF element
lambdaLEN = zeros(nelem,1) ; % Scale length paramters for each EIF element
D_Wecm = cell(nelem,1) ; % Diagonal matrix with weights
ndim = size(COOR,2);    

if ndim == 2
    D_AspMAT = cell(nelem,1) ; % Diagonal matrix   whose diagonal matrices are the "spin" operators AspMAT (A_sp)
elseif ndim == 3
    D_AspMAT = cell(nelem,ndim) ;
end



D_PdownsRBlROT = cell(nelem,1) ; % Diagonal matrix containing the matrices of each element relating increments/variations of the rotation vector
% in local coordinates with increments/variations of coarse-scale
% displacements 

% D_YcmpD = cell(nelem,1) ; % Diagonal matrix, residual for the compatibility condition
% D_Zfict = cell(nelem,1) ; % Diagonal matrix, linearization residual for the compatibility condition

%
% \XcALL \defeq  \colcuatro{\coldos{\XcALLbE{1}}{\zero}}{\coldos{\XcALLbE{2}}{\zero}}{\vdots}{\coldos{\XcALLbE{\nelemC}}{\zero}}
XcALL = cell(nelem,1) ;  % Global reference axes
XcALLloc = cell(nelem,1) ;  % Local reference axes 

WEIGHTSinteg.INTforces = cell(nelem,1) ;  % Determinant of the Jacobian matrix
WEIGHTSinteg.BodyForces = cell(nelem,1) ;  % Determinant of the J acobian matrix
%DATA =DefaultField(DATA,'IndexPermutationConnectivities',ones(size(MaterialType))) ;
DATA = DefaultField(DATA,'CriterionChooseParentDomain','max_Q_11') ;
DATA = DefaultField(DATA,'UNIFORM_SCALING_REFERENCE_ELEMENT',1) ;   % THIS IS THE DEFAULT OPTION, 11-March-2023
%   Criterion  for choosing the connectivity of the Element  (Cit might be changed within the code in order to minimize the dilatational component
% )
DATA.MESH.IndexECMpoints_per_element = cell(nelem,1) ;  % For instance, if we have two coarse-scale elements, with 6 and 10 integration points
% then IndexECMpoints_per_element{1} = 1:6,  and
% IndexECMpoints_per_element{2} = 7:16
% -----------------------------------------------------------------------
DATA.MESH.IndexDOFS_per_element = cell(nelem,1) ; % The same as above, but for DOFs


% CECM INTEGRATION RULE FOR NONLINEAR STRESSES (oPTION IMPLEMENTED ON
% 7-Jan-2024), see
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
if isfield(PROPMAT(1).EIFE_prop.INTforces,'CECM_ONLY_FOR_NONLINEAR_STRESSES')
    DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = PROPMAT(1).EIFE_prop.INTforces.CECM_ONLY_FOR_NONLINEAR_STRESSES ;
    disp('---------- CECM ONLY FOR NONLINEAR PART OF STRESSES -----')
    disp('Coarse-scale stiffness matrix   defined in  PROPMAT(imat).EIFE_prop.Kcoarse  (parent domain ) uses the same elastic constant' )
    disp('as thosed used in training... Do  not change them !  !!!! ')
else
    DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES = 0 ;
end

% THIS CO-ROTATIONAL APPROACH SHOULD BE HANDLED ANY NUMBER OF BUBBLE DOFS PER
% ELEMENT (IN CONTRAST TO B_N_matricesEIFEbub.m )
% ------------------------------------------------------------------------
BUBBLE_DOFS = zeros(size(PROPMAT)); % Number of BUBBLE_DOFs per type of element
for imat = 1:length(PROPMAT)
    BUBBLE_DOFS(imat)= length(PROPMAT(imat).EIFE_prop.INFO.DOFsBUB) ;
    PROPMAT(imat) = DefaultField(PROPMAT(imat),'CriterionChooseParentDomain','MAX_Q_11') ;
end
DATA = DefaultField(DATA,'DO_NOT_COMPUTE_NSHAPE',0) ; 
%DATA.DO_NOT_COMPUTE_NSHAPE  = 0 ; 
nDOFSbub = unique(BUBBLE_DOFS) ;
if length(nDOFSbub) ~= 1
    disp('This is a temporary amendment (28-Oct-2024). The N-shape matrix should be assembled using the same procedure as the B-matrix')
 %   error('The number of bubble modes is to be the same for all elements in this implementation')
    DATA.DO_NOT_COMPUTE_NSHAPE = 1; 
end


disp('-----------------------------------------------------------')
disp('Determination parent domains + Bmat, Nmat matrices')
disp('-----------------------------------------------------------')

if ~isempty(CN_SUPPORT)
    disp(['...This problem has uncoupled meshes...It is necessary to project onto a  support mesh for post-process']);
    nnodeE_support = size(CN_SUPPORT,2) ;   % Number of nodes support mesh (4, for linear quadrilateral)
    disp(['Number of nodes per element support mesh = ',num2str(nnodeE_support)]) ;
    nsupMESH = size(CN,2)/nnodeE_support ;  % Number of uncoupled meshes, up  to 26-Apr-2024, just 2
    disp(['Number of uncoupled meshes = ',num2str(nsupMESH)]) ;
else
    nnodeE_support = size(CN,2) ;
end

%DATA.PERMUT = PermutationConnectivities(TypeElement,size(CN,2)) ;
DATA.PERMUT = PermutationConnectivities(TypeElement,nnodeE_support) ;  % JAHO, 26-aPR-2024
% DATA.nnodeE_geometry may be different than size(CN,2)

% THESE ARE THE INDICES FOR ASSEMBLing in the online stage THE diagonal matrix
% containing the
% extended  ROTATION MATRICES of all the elements
% ----------------------------------------------------------------
number_accumulated_DOFs =0 ;
INDEXsparseROTmat.ROWS.DOFsB = cell(nelem,1) ;
INDEXsparseROTmat.ROWS.DOFsBUB = cell(nelem,1) ;
INDEXsparseROTmat.COLS.DOFsB = cell(nelem,1) ;
INDEXsparseROTmat.COLS.DOFsBUB = cell(nelem,1) ;
NumberBubbleDOFS_perELEMENT = zeros(nelem,1) ;

% THESE ARE THE INDICES FOR ASSEMBLing THE GLOBAL ROTATION MATRICES
% ----------------------------------------------------------------
INDEXsparseFINT.ROWS = cell(nelem,1) ;
INDEXsparseFINT.COLS = cell(nelem,1) ;



ACUM_ecm_points = 1;

for e = 1:nelem % Loop over coarse-scale elements
    disp(['e=',num2str(e)])
    CNloc = CN(e,1:nnodeE_support) ;   %  Nodes of element "e" (ordered by the mesher, in this case GID).
    
    EIFEoper_all = PROPMAT(MaterialType(e)).EIFE_prop ; % Properties EIF object ("element")
    DATA.CriterionChooseParentDomain =PROPMAT(MaterialType(e)).CriterionChooseParentDomain ;
    %  [CNnew,TRANSF_COORD{e},Vrot_all{e}] = ParentDomainSearchEIFE(EIFEoper_all,COOR,CNloc,TypeElement,DATA) ;
    [CNnew,TRANSF_COORD{e},PERMUT_chosen] = ParentDomainSearchEIFEcorot(EIFEoper_all,COOR,CNloc,TypeElement,DATA) ; % 26-Apr-2024
    
    
    QrotINI{e} = TRANSF_COORD{e}.QrotINI ; % INITIAL ROTATION MATRIX
    lambdaLEN(e) = TRANSF_COORD{e}.lambdaLEN ; % SCaling length parameters
    % Expanded rotation matrix (repeated for each coarse-scale node + identy for bubble DOFs)
    [~,iROWSloc,jCOLSloc] = QrotINIall_Element(EIFEoper_all,QrotINI{e}) ;
    
    % Global indices, assemblying diagonal rotation matrices
    INDEXsparseROTmat.ROWS.DOFsB{e} = iROWSloc.DOFsB + number_accumulated_DOFs;
    INDEXsparseROTmat.ROWS.DOFsBUB{e} = iROWSloc.DOFsBUB(:) +number_accumulated_DOFs;
    INDEXsparseROTmat.COLS.DOFsB{e} = jCOLSloc.DOFsB +number_accumulated_DOFs;
    INDEXsparseROTmat.COLS.DOFsBUB{e} = jCOLSloc.DOFsBUB(:) +number_accumulated_DOFs;
    NumberBubbleDOFS_perELEMENT(e) = length(iROWSloc.DOFsBUB ) ;
    
    nDOFS_element =    EIFEoper_all.INFO.DOFsBUB(end);
    DATA.MESH.IndexDOFS_per_element{e} = (number_accumulated_DOFs+1):(number_accumulated_DOFs+nDOFS_element);
    
    % Global indices, internal forces
    [ INDEXsparseFINT.ROWS{e},INDEXsparseFINT.COLS{e}] = FintCunassemb_Element(e,number_accumulated_DOFs,EIFEoper_all);
        number_accumulated_DOFs = number_accumulated_DOFs +  nDOFS_element  ;

   
    
    
    
    %   Rotated and scaled coordinates of  coarse-scale element
    %   \coldos{\XcALLbE{1}}{\zero} where
    %   \XcALLbE{e} \defeq \lambdaLENe{e} \DiagONEf{\QrotINIe{e}} \XcREFallBe{e}
    
    % MESH.XcREFallBe is computed in
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/EIFE_operatorsBUB.m
    
    XcALLbE =   lambdaLEN(e)*QrotINI{e}*(EIFEoper_all.MESH.XcREFallBe');
    XcALL{e} = [XcALLbE(:); zeros(NumberBubbleDOFS_perELEMENT(e),1)] ;
 
%  \XcALLloc \defeq  \colcuatro{\coldos{\XcALLbLOCe{1}}{\zero}}{\coldos{\XcALLbLOCe{2}}{\zero}}{\vdots}{\coldos{\XcALLbLOCe{\nelemC}}{\zero}}
%  \XcALLbE{e} = \DiagCb{\QrotINIe{e}}   \XcALLbLOCe{e}  
%   \XcALLbLOCe{e} =  \lambdaLENe{e}  \XcREFallBe{e}
    
    XcALLbLOCe =   lambdaLEN(e)*(EIFEoper_all.MESH.XcREFallBe');
    XcALLloc{e} = [XcALLbLOCe(:); zeros(NumberBubbleDOFS_perELEMENT(e),1)] ;
    
    
  % SPIN OPERATOR 
  [D_A,DOFsROT,DOFsTR] =  SpinOperator(ndim,size(XcALLbE,2),NumberBubbleDOFS_perELEMENT(e)) ; 
  
  if ndim == 2
      D_AspMAT{e} = D_A; 
  else
      D_AspMAT(e,:) = D_A ; 
  end
    
    % -----------------------------------------------------------------------
     % DOWNSCALING OPERATOR, RIGID BODY, ROTATIONAL COMPONENT, D_PdownsRBlROT
     % -------------------------------------------
    D_PdownsRBlROT{e} = D_PdownsRBlROT_COROT(EIFEoper_all,lambdaLEN(e),DOFsROT,DOFsTR) ;  
    
    
    
    
    if ~isempty(CN_SUPPORT)
        CN_SUPPORT(e,:) =   CN_SUPPORT(e,PERMUT_chosen) ;
        if nsupMESH == 2
            CN1 = CN(e,1:nnodeE_support)  ;  CN2 = CN(e,nnodeE_support+1:end) ;
            CN(e,:) = [CN1(PERMUT_chosen),CN2(PERMUT_chosen)] ;
        else
            error('Option not implemented')
        end
    else
        CN(e,:) = CNnew ;
    end
    % B-MATRIX, weights integration forces
    % --------------------------------------------------------
    [BmatIst_e,WEIGHTSinteg.INTforces{e},KcLINlocE] = Bmat_weights_EIFEbubCOROT(DATA,TRANSF_COORD{e},EIFEoper_all) ;
    
    nstrain = size(BmatIst_e,1)/length(WEIGHTSinteg.INTforces{e}) ;
    % Diagonal matrix
    D_Wecm{e} = repmat(WEIGHTSinteg.INTforces{e},1,nstrain) ;
    
    
    
    DATA.MESH.IndexECMpoints_per_element{e} = ACUM_ecm_points:(ACUM_ecm_points+length(WEIGHTSinteg.INTforces{e}) -1) ;
    ACUM_ecm_points  = ACUM_ecm_points + length(WEIGHTSinteg.INTforces{e}) ;
    
    
    D_Wecm{e} = D_Wecm{e}';
    D_Wecm{e} = D_Wecm{e}(:) ;
    D_Wecm{e} = diag(sparse(D_Wecm{e})) ;
    %-----------------------
    % Linear stiffness matrix
    D_KcLINloc{e} = sparse(KcLINlocE) ;
    D_BmatIst{e} = sparse(BmatIst_e) ;
    %  D_PdownsRB{e}  = sparse(EIFEoper_all.RECONSTRUCTION.RB_DISP.coeff) ;
    
    % Compatibility equation for rotation
    % Matrix appearing in the expression for the residual
    %   \DiagC{\YcmpD} \defeq \diagOL{\YcmpDe{1}}{\YcmpDe{2}}{\cdots}{\YcmpDe{\nelemC}}
    %
    %     D_YcmpD{e}  = sparse(EIFEoper_all.OPER.YcmpD) ;
    %     % For linearization of the compatibility residual
    %      D_Zfict{e}  = sparse(EIFEoper_all.OPER.Zfict) ;
    
    % N-matrix, weights integration body/inertial forces
    if DATA.DO_NOT_COMPUTE_NSHAPE == 0
    [Nmat{e},WEIGHTSinteg.BodyForces{e}] = Nmat_weights_EIFEbubCOROT(DATA,TRANSF_COORD{e},EIFEoper_all) ; 
    end
    
end


% ASSEMBLY EXTENDENDED ROTATION MATRICES (JUST INDICES )
% -------------------------------------------------------------------
QrotINI = cell2mat(QrotINI) ;
INDEXsparseROTmat.ROWS.DOFsB = cell2mat(INDEXsparseROTmat.ROWS.DOFsB) ;
INDEXsparseROTmat.COLS.DOFsB = cell2mat(INDEXsparseROTmat.COLS.DOFsB) ;
INDEXsparseROTmat.ROWS.DOFsBUB = cell2mat(INDEXsparseROTmat.ROWS.DOFsBUB) ;
INDEXsparseROTmat.COLS.DOFsBUB = cell2mat(INDEXsparseROTmat.COLS.DOFsBUB) ;
INDEXsparseROTmat.nelem = nelem;
CHECK_DIAG_MATRIX = 0 ;
if CHECK_DIAG_MATRIX == 1
    D_QrotINIall = QrotINIall_GLOBAL(QrotINI,INDEXsparseROTmat) ;
end

% ASSEMBLY MATRIX ROTATION MATRIX STORED BY ROWS (ONLY 3D)
if ndim ==3 
    [irows,icols] = IndexesQrotROWS(nelem,ndim) ; 
    INDEXsparseROTelem_byROWS.IROWS = irows ;  
    INDEXsparseROTelem_byROWS.ICOLS = icols ;  
    CHECK_DIAG_MATRIX = 0 ;  
    if CHECK_DIAG_MATRIX == 1
    D_QrotINI_byROWS = Qrot_by_rows_element(QrotINI,INDEXsparseROTelem_byROWS) ;
    end
else
    INDEXsparseROTelem_byROWS = [] ; 
end




% INTERNAL FORCES
INDEXsparseFINT.ROWS = cell2mat(INDEXsparseFINT.ROWS);
INDEXsparseFINT.COLS = cell2mat(INDEXsparseFINT.COLS);

%
XcALL = cell2mat(XcALL) ;
XcALLloc = cell2mat(XcALLloc) ;

% Converting BmatIst into a block diagonal matrix
D_BmatIst = blkdiag(D_BmatIst{:}) ;
% Matrix weights for internal forces (nonlinear part)
D_Wecm = blkdiag(D_Wecm{:}) ;
% Linear stiffness matrix
D_KcLINloc= blkdiag(D_KcLINloc{:}) ;
D_PdownsRBlROT = blkdiag(D_PdownsRBlROT{:}) ;

if ndim ==2
D_AspMAT = blkdiag(D_AspMAT{:}) ;
elseif ndim == 3
    D_AspMAT_glo = cell(1,ndim) ; 
    for idim = 1:ndim 
        D_AspMAT_loc = D_AspMAT(:,idim) ; 
        D_AspMAT_glo{idim} = blkdiag(D_AspMAT_loc{:});         
    end
    D_AspMAT = D_AspMAT_glo ; 
else
    error('Option not implemented')
end


% Boolean operator (for assembly)
[LboolCall,~] =   BooleanOperatorCoarse(COOR,CN,number_accumulated_DOFs,NumberBubbleDOFS_perELEMENT) ;

 
%
%
% GHOST DOFS AND EXTENDED VARIABLES
MESHextended = [] ;
[MESHextended.NDOFS_pernode,MESHextended.COOR,MESHextended.CN,MESHextended.DOFS_TO_KEEP,...
    MESHextended.DOFS_bLOC,MESHextended.DOFS_bubLOC,MESHextended.DOFS_ghost] = GhostDOFs(COOR,CN,DATA,nDOFSbub) ;
%
DATA.MESHextended = MESHextended ;

% maxDOF = max(bbb,size(MESHextended.COOR,2));
% % Maximum number of DOFs per node
% disp(['----------------------------------------------'])
% disp(['Maximum Number of DOFs per node'])
% disp(['DOFmax = ',num2str(maxDOF)])
% disp(['----------------------------------------------'])
% DATA.maxNDOFperNODE = maxDOF ;




