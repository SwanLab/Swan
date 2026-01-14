function [D_BmatIst,Nmat,WEIGHTSinteg,TRANSF_COORD,CN,DATA,D_KcLINloc,CN_SUPPORT,QrotINI,lambdaLEN,D_Wecm,INDEXsparseROTmat,...
    LboolCall,XcALL,D_PdownsRB] = B_N_matricesEIFEbubCOROT(COOR,CN, PROPMAT,MaterialType,DATA,TypeElement,...
    CN_SUPPORT)
% EIFE METHOD, DETERMINATION OF B and N matrices  --------------
% Adaptation of B_N_matricesEIFEbub.m
% CO-ROTATIONAL APPROACH, 22-OCT-2O24, BALMES 185, BARCELONA
% --------------------------------------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
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
D_PdownsRB = cell(nelem,1) ; % Diagonal matrix, downscaling operator for rigid body amplitudes
% 
% \XcALL \defeq  \colcuatro{\coldos{\XcALLbE{1}}{\zero}}{\coldos{\XcALLbE{2}}{\zero}}{\vdots}{\coldos{\XcALLbE{\nelemC}}{\zero}}
XcALL = cell(nelem,1) ; 

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

nDOFSbub = unique(BUBBLE_DOFS) ;
if length(nDOFSbub) ~= 1
    disp('This is a temporary amendment (28-Oct-2024). The N-shape matrix should be assembled using the same procedure as the B-matrix')
    error('The number of bubble modes is to be the same for all elements in this implementation')
    
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


number_accumulated_DOFs =0 ;
INDEXsparseROTmat.ROWS.DOFsB = cell(nelem,1) ;
INDEXsparseROTmat.ROWS.DOFsBUB = cell(nelem,1) ;
INDEXsparseROTmat.COLS.DOFsB = cell(nelem,1) ;
INDEXsparseROTmat.COLS.DOFsBUB = cell(nelem,1) ;
NumberBubbleDOFS_perELEMENT = zeros(nelem,1) ; 
ACUM_ecm_points = 1; 
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
    INDEXsparseROTmat.ROWS.DOFsBUB{e} = iROWSloc.DOFsBUB +number_accumulated_DOFs;
    INDEXsparseROTmat.COLS.DOFsB{e} = jCOLSloc.DOFsB +number_accumulated_DOFs;
    INDEXsparseROTmat.COLS.DOFsBUB{e} = jCOLSloc.DOFsBUB +number_accumulated_DOFs;
    NumberBubbleDOFS_perELEMENT(e) = length(iROWSloc.DOFsBUB ) ; 
    
   nDOFS_element =    EIFEoper_all.INFO.DOFsBUB(end); 
    
   DATA.MESH.IndexDOFS_per_element{e} = (number_accumulated_DOFs+1):(number_accumulated_DOFs+nDOFS_element);
    
    number_accumulated_DOFs = number_accumulated_DOFs +  nDOFS_element  ;
    
   
    
    %   Rotated and scaled coordinates of  coarse-scale element 
    %   \coldos{\XcALLbE{1}}{\zero} where
    %   \XcALLbE{e} \defeq \lambdaLENe{e} \DiagONEf{\QrotINIe{e}} \XcREFallBe{e}
    
    
    % MESH.XcREFallBe is computed in 
   % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/EIFE_operatorsBUB.m
    
    
     XcALLbE =   lambdaLEN(e)*QrotINI{e}*(EIFEoper_all.MESH.XcREFallBe');
    XcALL{e} = [XcALLbE(:); zeros(NumberBubbleDOFS_perELEMENT(e),1)] ; 
    
    
    
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
    D_PdownsRB{e}  = sparse(EIFEoper_all.RECONSTRUCTION.RB_DISP.coeff) ;
     
    % N-matrix, weights integration body/inertial forces
    [Nmat{e},WEIGHTSinteg.BodyForces{e}] = Nmat_weights_EIFEbubCOROT(DATA,TRANSF_COORD{e},EIFEoper_all) ;

    
    
    
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

%
XcALL = cell2mat(XcALL) ; 
% Converting BmatIst into a block diagonal matrix 
D_BmatIst = blkdiag(D_BmatIst{:}) ; 
% Matrix weights for internal forces (nonlinear part)
D_Wecm = blkdiag(D_Wecm{:}) ; 
 % Linear stiffness matrix     
D_KcLINloc= blkdiag(D_KcLINloc{:}) ; 
% Boolean operator (for assembly)
[LboolCall,XcALL_absol] =   BooleanOperatorCoarse(COOR,CN,number_accumulated_DOFs,NumberBubbleDOFS_perELEMENT) ; 

% disp('borrar esto....')
% XcALL =  XcALL_absol ; 

% Downscaling operator, rigid body 
D_PdownsRB= blkdiag(D_PdownsRB{:}) ; 



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




