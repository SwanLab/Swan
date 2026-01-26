function [EIFEoper,OPERFE,MESHdom,GEOMATRICES,MODES] = EIFE_operatorsFUN(INPUT_PROBLEMS,DATAoffline)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% Empirical Interscale Finite Element Method (EIFEM)
% ------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end

DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',[]) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'INTERNAL_FORCES',0) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'MASS_MATRIX',0) ;

DATAcommon  = feval(INPUT_PROBLEMS) ;
NAME_BASE =[DATAcommon.NameParamStudyLOC,'_param_'];
NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];

% load offline information
load(NAMEOFFLINEstore,'MODES','GEOMATRICES','MESHdom','DATA_misc','MATPRO','OPERFE','CECM','DATAcommon')  ;


% 1) Hdef
%     \begin{equation}
% \label{eq:deHs1}
%  \Hdef{e}{} \defeq  \PsiDEF{e}{}^T \PhiDEF{e}{}
f = MESHdom.faceDOFSall ;
% Work done by the self-equilibrated modes over deformational displacements
%
Hdef = MODES.PsiDEFf'*MODES.PhiDEF(f,:) ;

% 2) Hrb
% ----------
% \begin{equation}
% \label{eq:deHs1bis}
%  \Hrb{e}{} \defeq  \PsiRB{e}{}^T \PhiDEF{e}{}

% Work done by the resultant modes over deformational displacements
Hrb = MODES.PsiRBf'*MODES.PhiDEF(f,:) ;


% 3) Grb
% \begin{equation}
% \label{eq:Gdefinition}
%  \Grb{e}{} \defeq \PsiRB{e}{}^T \PhiRB{e}{} =    \PsiRB{e}{\faceDOFS{}{1}}^T \PhiRB{e}{\faceDOFS{}{1}} +   \PsiRB{e}{\faceDOFS{}{2}}^T \PhiRB{e}{\faceDOFS{}{2}}    , \hspace{1cm}       e= 1,2 \ldots \nDOM
% \end{equation}

% Work done by the resultant modes over rigid-body modes

Grb = MODES.PsiRBf'*MODES.PhiRB(f,:) ;

% 4) Work done by self-equilibrated forces over interface displacements
%  \Tdef{e}{} \defeq \PsiDEFf{e}{}^T \Dintf{e}{} \Vall{e}{}

Tdef = MODES.PsiDEFf'*MODES.Vall ;

% Deformational part
coeff = (MODES.PhiRB(f,:)'*GEOMATRICES.Mintf*MODES.PhiRB(f,:))\(MODES.PhiRB(f,:)'*GEOMATRICES.Mintf*MODES.Vall);
VallDEF = MODES.Vall - MODES.PhiRB(f,:)*coeff ;
DATALOC.TOL  = 1e-4;
[VallDEF,SS,VV] = WSVDT(VallDEF, (GEOMATRICES.Mintf),DATALOC) ;
CONDITION_WELLPOSEDNESS =  MODES.PsiDEFf'*VallDEF ;
[UUU,SSS,VVV]  = SVDT(CONDITION_WELLPOSEDNESS);
disp('Singular values WORK DONE AT THE INTERFACES')
disp(num2str(SSS'))

% 4) Work done by resultant forces over interface displacements

Trb = MODES.PsiRBf'*MODES.Vall ;

% Downscaling operators
% ---------------------
%  \PdownsDEF{e}{}  =  \Hdef{e^{-1}}{} \Tdef{e}{},

PdownsDEF  = Hdef\Tdef ;
%
%   \begin{equation}
%   \PdownsRB{e}{}  =   \Grb{e}{}^{-1} \Par{     \Trb{e}{} -\Hrb{e}{}  \Hdef{e^{-1}}{} \Tdef{e}{}  },   \hspace{1cm} e= 1,2 \ldots \nDOM
%  \end{equation}

PdownsRB  = Grb\(Trb - Hrb*(Hdef\Tdef) )  ;

if isfield(MATPRO,'Celas')
    MATPRO.celasglo = MATPRO.Celas ;
end

EIFEoper = [] ;

if ~isempty(MATPRO.celasglo)
    %     disp('Checking accuracy CECM points ')
    %     disp('-------------------------------------------')
    %   disp('Coarse-scale stiffness matrix ')
    celastST = MATPRO.celasglo ;
    nF = size(MATPRO.celasglo,2) ;
    for icomp = 1:nF
        icol = icomp:nF:size(celastST,1) ;
        celastST(icol,:) = bsxfun(@times,celastST(icol,:),OPERFE.wSTs) ;
    end
    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
    Kfine = OPERFE.Bst'*(celastST*OPERFE.Bst);
    
    
    
    
        % CHECKING ORTHOGONALITY CONDITIONS
    % ----------------------------------
    CHECK_ORTH = 1;
    if CHECK_ORTH == 1
        disp('This is just for checking purposes')
        
        disp('F = Kfine*PhiDEF')
        Fdef = Kfine*MODES.PhiDEF ; 
        ALLDOFS = 1:size(Fdef,1) ; 
        INT_DOFS = setdiff(ALLDOFS,f) ; 
        perc_inter = norm(Fdef(INT_DOFS,:),'fro')/ norm(Fdef,'fro')*100; 
        disp(['Relat. Norm F_def(INTERIOR,:) = ',num2str(perc_inter),' %'])
        
        
        
        
    end
    
    Kred = MODES.PhiDEF'*Kfine*MODES.PhiDEF;
    
    Kcoarse= PdownsDEF'*Kred*PdownsDEF ;
    
    %     dbstop('98')
    %         disp('Borrar esto (temporal). Guarda información')
    %     PhiDEF=  MODES.PhiDEF ;
    %     save('DATA_INFLUENCE_size/DATA_p4.mat','Kcoarse','PhiDEF')
    
    
    % Mass matrix
    
    Mred_rb  = MODES.PhiRB'*OPERFE.M*MODES.PhiRB ;
    
    Mcoarse_rb = PdownsRB'*Mred_rb*PdownsRB ;
    
    Mred_def  = MODES.PhiDEF'*OPERFE.M*MODES.PhiDEF ;
    
    Mcoarse_def = PdownsDEF'*Mred_def*PdownsDEF ;
    
    Mcoarse = Mcoarse_rb + Mcoarse_def ;
    
    
    save('Kcoarse_homog.mat','Kcoarse','Mcoarse')
    
    
    EIFEoper.Kcoarse = Kcoarse ; 
    EIFEoper.Mcoarse = Mcoarse ; 
    
end

% ------------------------------------------------------------------------------------
% B-matrices, weights and points for integrating internal forces  (Right-hand side)
% -------------------------------------------------------------------------------------
EIFEoper.INTforces = [] ;
%
% METHOD =1;
%
% if METHOD == 1
% Version before March 2nd, 2023, 10:48, see comments about the change in
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/04_DistortedElement.mlx
%     [Bmat,wECM,xECM,MATPROpoints,RECONSTPK2,index_elements,BmatRED] ...
%         = BmatricesEIFE_v1(MODES.PhiDEF,PdownsDEF,CECM.INTERNAL_FORCES,OPERFE.Bst,...
%         DATAoffline,DATA_misc,MATPRO,MODES.BasisStwo,MESHdom) ;

[Bmat,wECM,xECM,MATPROpoints,RECONSTPK2,index_elements,BmatRED,Bgrad,PhiDEFelem] ...
    = BmatricesEIFE_v2(MODES.PhiDEF,PdownsDEF,CECM.INTERNAL_FORCES,OPERFE.Bst,...
    DATAoffline,DATA_misc,MATPRO,MODES.BasisStwo,MESHdom) ;


EIFEoper.INTforces.BmatRED = BmatRED;
EIFEoper.INTforces.Bgrad_cell = Bgrad;
EIFEoper.INTforces.PhiDEFelem_cell = PhiDEFelem;

EIFEoper.INTforces.weights = wECM;
EIFEoper.INTforces.posgp = xECM;
EIFEoper.INTforces.MATPRO = MATPROpoints;
EIFEoper.INTforces.MaterialType = MESHdom.MaterialType(index_elements);

%else
%     % Before 1-March-2023
%     % Separated representation of the B-matrix
%     [Bmat,wECM,xECM,MATPROpoints,RECONSTPK2,index_elements,Bgrad,PhiDEFelem] ...
%         = BmatricesEIFE(MODES.PhiDEF,PdownsDEF,CECM.INTERNAL_FORCES,OPERFE.Bst,DATAoffline,DATA_misc,MATPRO,MODES.BasisStwo,MESHdom) ;
%
%     EIFEoper.INTforces.Bmat = Bmat;
%     EIFEoper.INTforces.Bgrad_cell = Bgrad;
%     EIFEoper.INTforces.PhiDEFelem_cell = PhiDEFelem;
%
%     EIFEoper.INTforces.weights = wECM;
%     EIFEoper.INTforces.posgp = xECM;
%     EIFEoper.INTforces.MATPRO = MATPROpoints;
%     EIFEoper.INTforces.MaterialType = MESHdom.MaterialType(index_elements);
%
% end




% CHECK ACCURACY APPROXIMATING STIFFNESS MATRIX
EIFEoper.MESH.nstrain         = DATA_misc.MESH.nstrain ;
DATA_misc  = DefaultField(DATA_misc,'FESHAPE_coarse_elem_transf_coord',[]) ;

DATA_misc.FESHAPE_coarse_elem_transf_coord  = DefaultField(DATA_misc.FESHAPE_coarse_elem_transf_coord,'ShapeFunction','') ;

switch  DATA_misc.FESHAPE_coarse_elem_transf_coord.ShapeFunction
    case 'ShapeFunctionFE'
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/FictInterfaceDISP_QUADlin.m
        DATAlocSHAPE = DATA_misc.FESHAPE_coarse_elem_transf_coord ;
        %     % Gauss points body forces
        %     [Nshape_transf, Bgrad_transf ]=    ShapeFunctionFE([],EIFEoper.BodyForces.posgp,DATAlocSHAPE.DATA_ShapeFunctionFE) ;
        %     EIFEoper.BodyForces.Nshape_transf = Nshape_transf  ;
        %     EIFEoper.BodyForces.Bgrad_transf = Bgrad_transf  ;
        
        % Gauss points internal forces
        if ~isempty(DATAlocSHAPE.DATA_ShapeFunctionFE.DATAshape)
            [Nshape_transf, Bgrad_transf ]=    ShapeFunctionFE([],EIFEoper.INTforces.posgp,DATAlocSHAPE.DATA_ShapeFunctionFE) ;
            EIFEoper.INTforces.Nshape_transf = Nshape_transf  ;
            EIFEoper.INTforces.Bgrad_transf = Bgrad_transf  ;
        end
        
    case 'ShapeFunctionFEtri'
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/FictInterfaceDISP_QUADlin.m
        DATAlocSHAPE = DATA_misc.FESHAPE_coarse_elem_transf_coord ;
        %     % Gauss points body forces
        %     [Nshape_transf, Bgrad_transf ]=    ShapeFunctionFE([],EIFEoper.BodyForces.posgp,DATAlocSHAPE.DATA_ShapeFunctionFE) ;
        %     EIFEoper.BodyForces.Nshape_transf = Nshape_transf  ;
        %     EIFEoper.BodyForces.Bgrad_transf = Bgrad_transf  ;
        
        % Gauss points internal forces
        [Nshape_transf, Bgrad_transf ]=    ShapeFunctionFEtri([],EIFEoper.INTforces.posgp,DATAlocSHAPE.DATA_ShapeFunctionFE) ;
        EIFEoper.INTforces.Nshape_transf = Nshape_transf  ;
        EIFEoper.INTforces.Bgrad_transf = Bgrad_transf  ;
    otherwise
        EIFEoper.INTforces.Nshape_transf = [] ;
        EIFEoper.INTforces.Bgrad_transf = [] ;
        
end

if ~isempty(MATPRO.celasglo)
    %  disp('Coarse-scale stiffness matrix APPROXIMATED BY ECM POINTS/WEIGHTS ')
    
    
    
    celastST = MATPROpoints.celasglo ;
    nF = size(MATPROpoints.celasglo,2) ;
    for icomp = 1:nF
        icol = icomp:nF:size(celastST,1) ;
        celastST(icol,:) = bsxfun(@times,celastST(icol,:),wECM) ;
    end
    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
    KCOARSE_ECM = Bmat'*(celastST*Bmat);
    
    ERROR_approx = norm(KCOARSE_ECM-Kcoarse,'fro')/norm(Kcoarse,'fro') ;
    disp(['Error approximation stiffness with respect integration full FE points: ',num2str(ERROR_approx)])
    
    
    %     % COMPARISON WITH ALTERNATIVE METHOD (Separated contribution)
    %     Xe = DATA_misc.FESHAPE_coarse_elem_transf_coord.COOR' ;
    %     Ke = ComputeKeMatrix_multiNV(EIFEoper,Xe) ;
    %
    %        ERROR_approx = norm(KCOARSE_ECM-Ke,'fro')/norm(KCOARSE_ECM,'fro') ;
    %     disp(['Error approximation stiffness with respect to vectorized approach: ',num2str(ERROR_approx)])
end

% ------------------------------------------------------------------------------------
% N-matrices, weights and points for integrating body external forces  (Right-hand side)
% -------------------------------------------------------------------------------------
EIFEoper.BodyForces = [] ;

% switch DATAcommon.TypeFunctionDisplacementInterfaces
%     case 'PLATE_QUAD_LINEAR'
%         DATA_misc.MESH.ndim =5 ;
% end

[Nmat,wECM,xECM,MATPROpoints,index_elements] ...
    = NmatBODY_EIFE(MODES.PhiDEF,MODES.PhiRB,PdownsDEF,PdownsRB,CECM.MASS_MATRIX,OPERFE.Nst,DATAoffline,DATA_misc,MATPRO) ;
EIFEoper.BodyForces.Nmat = Nmat;
% EIFEoper.BodyForces.NmatDEFred = Nmat_def_red;  % Deformational part
% EIFEoper.BodyForces.NmatRBred  = Nmat_rb_red;    % Rigid body part

EIFEoper.BodyForces.weights = wECM;
EIFEoper.BodyForces.posgp = xECM;
EIFEoper.BodyForces.MATPRO = MATPROpoints;
EIFEoper.BodyForces.MaterialType = MESHdom.MaterialType(index_elements);


% CHECK ACCURACY APPROXIMATING mass MATRIX

if ~isempty(MATPRO.celasglo)
    % disp('Coarse-scale mass matrix APPROXIMATED BY ECM POINTS/WEIGHTS ')
    % 5.  MASS MATRIX
    % -----------------------
    densGLO = repmat(MATPROpoints.dens',length(MATPROpoints),1) ;
    densGLO = densGLO(:) ;
    densGLO_W = CompWeightDiag(densGLO.*wECM,DATA_misc.MESH.ndim) ;
    NstW = densGLO_W*Nmat ;
    Mcoarse_ECM = NstW'*Nmat ;
    
    ERROR_approx = norm(Mcoarse-Mcoarse_ECM,'fro')/norm(Mcoarse,'fro') ;
    disp(['Error approximation mass matrix with respect integration full FE points: ',num2str(ERROR_approx)])
end


%% RECONSTRUCTION OPERATORS
% ------------------------------
% Format =  BASIS*coeff

% RIGID-BODY DISPLACEMENTS
EIFEoper.RECONSTRUCTION.RB_DISP.coeff = PdownsRB ;
EIFEoper.RECONSTRUCTION.RB_DISP.BASIS = MODES.PhiRB ;
% DEFORMATIONAL DISPLACEMENTS
EIFEoper.RECONSTRUCTION.DEF_DISP.coeff = PdownsDEF ;
EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS = MODES.PhiDEF ;
% STRAINS
EIFEoper.RECONSTRUCTION.STRAINS.coeff = PdownsDEF ;
EIFEoper.RECONSTRUCTION.STRAINS.BASIS = OPERFE.Bst*MODES.PhiDEF ;
EIFEoper.MODES.PsiDEFf = MODES.PsiDEFf ;
EIFEoper.MODES.PsiRBf = MODES.PsiRBf ;
EIFEoper.MODES.PhiDEF = MODES.PhiDEF ;
EIFEoper.MODES.PhiRB = MODES.PhiRB ;
EIFEoper.MODES.Vall = MODES.Vall ;
%EIFEoper.MODES.Bst_FE = OPERFE.Bst; 

% For stiffness matrix/internal forces
% -------------------------------------
EIFEoper.OPER.HdefINV_PsiDEFfT = Hdef\MODES.PsiDEFf' ;
% For mass matrix
Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT = Grb\(MODES.PsiRBf' - Hrb*(Hdef\MODES.PsiDEFf') )  ;
EIFEoper.OPER.Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT = Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT;


if ~isempty(MATPRO.celasglo)
    
    % RECONSTRUCTION BASIS MATRIX INCLUDE FOR STRESSES INCLUDE THE TERM
    % Celas*Bred*Ma
    
    celastST = ConvertBlockDiag(MATPRO.celasglo) ;
    
    RECONSTPK2.Celas_Bmat_PhiDEF =    celastST*EIFEoper.RECONSTRUCTION.STRAINS.BASIS;
    
    
    
end


%  PK2-STRESSES
EIFEoper.RECONSTRUCTION.PK2STRESS = RECONSTPK2 ;

% MESH PROPERTIES
% -----------------------
EIFEoper.MESH.COOR         = MESHdom.COOR ;
EIFEoper.MESH.CN           = MESHdom.CN ;
EIFEoper.MESH.TypeElement  = MESHdom.TypeElement ;
EIFEoper.MESH.MaterialType = MESHdom.MaterialType ;
EIFEoper.MESH.posgp         = DATA_misc.MESH.posgp ;

% ADDITIONAL INFORMATION
% ----------------------
EIFEoper.INFO.SMALL_STRAIN_KINEMATICS =  DATA_misc.SMALL_STRAIN_KINEMATICS ; % Small strain kinematics
EIFEoper.INFO.FILE_MATERIAL_DATA =  DATAcommon.FILE_MATERIAL_DATA ; % Material file
EIFEoper.INFO.FILE_MATERIAL_DATA =  DATAcommon.FILE_MATERIAL_DATA ; % Material file
EIFEoper.INFO.FESHAPE_coarse_elem_transf_coord =  DATA_misc.FESHAPE_coarse_elem_transf_coord ; % Material file

% Coordinates transformation using the FE shapefunctions
switch  DATA_misc.FESHAPE_coarse_elem_transf_coord.ShapeFunction
    case 'ShapeFunctionFE'
        % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/MultiHROM/FictInterfaceDISP_QUADlin.m
        %   DATAlocSHAPE = DATA_misc.FESHAPE_coarse_elem_transf_coord ;
        % Gauss points body forces
        if ~isempty(DATAlocSHAPE.DATA_ShapeFunctionFE.DATAshape)
            [Nshape_transf, Bgrad_transf ]=    ShapeFunctionFE([],EIFEoper.BodyForces.posgp,DATAlocSHAPE.DATA_ShapeFunctionFE) ;
            EIFEoper.BodyForces.Nshape_transf = Nshape_transf  ;
            EIFEoper.BodyForces.Bgrad_transf = Bgrad_transf  ;
        end
        %
        %      % Gauss points internal forces
        %     [Nshape_transf, Bgrad_transf ]=    ShapeFunctionFE([],EIFEoper.INTforces.posgp,DATAlocSHAPE.DATA_ShapeFunctionFE) ;
        %     EIFEoper.INTforces.Nshape_transf = Nshape_transf  ;
        %     EIFEoper.INTforces.Bgrad_transf = Bgrad_transf  ;
    case 'ShapeFunctionFEtri'
        [Nshape_transf, Bgrad_transf ]=    ShapeFunctionFEtri([],EIFEoper.BodyForces.posgp,DATAlocSHAPE.DATA_ShapeFunctionFE) ;
        EIFEoper.BodyForces.Nshape_transf = Nshape_transf  ;
        EIFEoper.BodyForces.Bgrad_transf = Bgrad_transf  ;
    otherwise
        EIFEoper.BodyForces.Nshape_transf = [] ;
        EIFEoper.BodyForces.Bgrad_transf = [] ;
        
end

% OPERATORS LATERAL SURFACES
EIFEoper.LateralM_operator = cell(size(MESHdom.BNDlateral)) ;
if ~isempty(MESHdom.BNDlateral)
    
    for ilateral = 1:length(MESHdom.BNDlateral)
        PROPlat = MESHdom.BNDlateral(ilateral) ;
        DOFsLOC = small2large(PROPlat.NODES,DATA_misc.MESH.ndim) ;
        % Linear quadrilateral shape function
        COORmid = DATA_misc.FESHAPE_coarse_elem_transf_coord.COOR(:,1:2) ;
        ORDER_POLYNOMIALS = ones(1,2) ;
        DATAshape = ShapeFunCoefficients(COORmid(:,1:2),ORDER_POLYNOMIALS) ;
        DATAlocSHAPE.DATAshape  = DATAshape;
        DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
        [Nshape, ~,~ ]=    ShapeFunctionFE([],PROPlat.COORrelA_global(:,1:2),DATAlocSHAPE) ;
        ndim = DATA_misc.MESH.ndim;
        Ne = sparse(ndim*size(Nshape,1),ndim*size(Nshape,2)) ;
        for idim = 1:ndim
            Ne(idim:ndim:end,idim:ndim:end) = Nshape ;
        end
        
        %         ndim = DATA_misc.MESH.ndim;
        Mdof = sparse(ndim*size(PROPlat.GeometricMassMatrix,1),ndim*size(PROPlat.GeometricMassMatrix,2)) ;
        for idim = 1:ndim
            Mdof(idim:ndim:end,idim:ndim:end) = PROPlat.GeometricMassMatrix ;
        end
        
        % F_nodes =Mdof*f_nodes  --> f_nodes = Distributed forces at the
        % fine-scale nodes of the lateral forces
        % F_nodes = Equivalent forces at the same nodes
        
        Mleft_RB = (MODES.PhiRB(DOFsLOC,:)*PdownsRB )'*Mdof; % RIGID BODY
        Mleft_DEF = (MODES.PhiDEF(DOFsLOC,:)*PdownsDEF )'*Mdof; % DEFORMATIONAL
        Mleft = Mleft_RB + Mleft_DEF;
        % Here    F_coarse_nodes = Mleft*f_nodes. F_coarse_nodes is the vector of coares-scale external forces (20 entries)
        % Finally
        EIFEoper.LateralM_operator{ilateral} = Mleft*Ne ;
        
        
        %
        
        
        %[nnodeBND, ndim] = size(MESHcoar_COOR) ;
        % Order of polynomial
        %
        %DATAshape = ShapeFunCoefficients(COORmid(:,1:2),ORDER_POLYNOMIALS) ;
        %DATAlocSHAPE.DATAshape  = DATAshape;
        %xLIM = [] ;
        %DATAlocSHAPE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
        %[Nshape, ~,~ ]=    ShapeFunctionFE(xLIM,COORbnd,DATAlocSHAPE) ;
        
        
        
        %         PROPlat = MESHdom.BNDlateral(ilateral) ;  % Properties surface
        %         DOFsLOC = small2large(PROPlat.NODES,DATA_misc.MESH.ndim) ;
        
        %
        %
        %
        %     Mcoarse_def = PdownsDEF'*Mred_def*PdownsDEF ;
        %
        %     Mcoarse = Mcoarse_rb + Mcoarse_def ;
        
    end
    
end



%
%  Xe = DATA_misc.FESHAPE_coarse_elem_transf_coord.COOR' ;
% Me = ComputeMeMatrixR_multi_test(Xe,EIFEoper,MODES.Vall) ;
%       % COMPARISON WITH ALTERNATIVE METHOD (Separated contribution)
%
%            ERROR_approx = norm(Mcoarse_ECM-Me,'fro')/norm(Mcoarse_ECM,'fro') ;
%          disp(['Error approximation mass matrix with respect to vectorized approach: ',num2str(ERROR_approx)])
%
