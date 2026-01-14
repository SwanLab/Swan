function  [PsiSEf,PhiDEF,GammaBUB] =   SelfEquilibratedModes_plas(PsiRBf,Mintf,AsnapREAC,DATALOC,DATAoffline,MESH,BasisUdeform,Mdom)
% COMPUTATION OF SELF-EQUILIBRATED MODES
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/04_GeneralTheory.mlx

if nargin == 0
    load('tmp1.mat')
end


% TEST: Dimension space spanned by the reactive forces (all snapshots).
% Just for testing purposes

% RELTOL = [0] ;
% DATAlocS = [] ;
% [U,S,V] = SRSVD(AsnapREAC,RELTOL,DATAlocS) ;

% No, they are not
% Basic/Additional tests


AsnapREACse = cell(size(AsnapREAC)) ;

b = MESH.faceDOFSall ;  % Boundary DOFs

% Self-equilibrated part (in principle, not necessary if one trains without external forces, just prescribed displacements)

for iproj = 1:length(AsnapREAC)
    [AsnapREACse{iproj},MintfINV] =  HprojSEf_operator(PsiRBf,Mintf,AsnapREAC{iproj}(b,:)) ;
end

% -----------------------------------------------------------------------
SNAPbasic = cell2mat(AsnapREACse(DATAoffline.INDEX_BASIC_TESTS)) ;
SNAPcompl = AsnapREACse(DATAoffline.INDEX_COMPL_TESTS) ;
% ---------------------------------------------------------


CHECK_BASIC_COMPLEMENTARY_MODES = 0 ;

if CHECK_BASIC_COMPLEMENTARY_MODES == 1
    disp('Determining  basic self-equilibrated modes')
    DATALOCww.TOL = 1e-8;
    %DATALOCww.Mchol = Mintfinv_chol ;
    [ PsISEfBS,Sbs,~,Mintfinv_chol] = WSVDT( SNAPbasic,MintfINV,DATALOCww) ;
    
    DATALOC.LEGEND_MODES_SE = 'Basic' ;
    PlotModesSE_SubdomainLevel(DATALOC,PsiRBf,PsISEfBS,MESH) ;
    
    
    
    %%%% METHOD FOR DETERMINING REMAINING MODES (ALL AT ONCE)
    disp('----------------------------------------------------------------------')
    disp(['Computing  self-equilibrated  modes via SVD (same tolerance as displacements)'])
    disp('----------------------------------------------------------------------')
    
    nproj = 1 + length(SNAPcompl) ;
    A  = cell(1,nproj) ;
    A{1} = Mintfinv_chol*SNAPbasic ;
    for iproj = 1:length(SNAPcompl)
        A{iproj+1} =  Mintfinv_chol*SNAPcompl{iproj} ;
    end
    TOL = DATAoffline.TOLSVD_complementary_modes*ones(size(SNAPcompl)) ;
    RELTOL = [1e-10,TOL] ;
    DATAlocS = [] ;
    [U,S,V] = SRSVD(A,RELTOL,DATAlocS) ;
    
    PsiSEf_viaSVD = Mintfinv_chol\U ;  % All modes (including basic ones)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
disp('----------------------------------------------------------------------')
disp(['Computing self-equilibrated modes via right-singular vectors displacements'])
disp('----------------------------------------------------------------------')
% A = U*S*V^T -->   U = (A*V)*inv(S)

disp(['This operation may be optimized...'])
PsiSEf_viaCOEFF = cell2mat(AsnapREACse)*BasisUdeform.V_upsilonDEF ;
% Normalization
DATALOCcc.Mchol = [] ;
[ Y,~,~,~] = WSVDT( PsiSEf_viaCOEFF,MintfINV,DATALOCcc) ;

if CHECK_BASIC_COMPLEMENTARY_MODES == 1
    
    % WHAT IS THE INTERSECTION BETWEEN PsiSEf_viaCOEFF AND PsiSEf_viaSVD
    disp(['Cosines angles between SVD-truncated space and the coefficient space'])
    disp(['We shall use  the COEFFICIENT SPACE '])
    [AA,SS,VV] = SVDT(PsiSEf_viaCOEFF) ;
    [BB,SS,VV] = SVDT(PsiSEf_viaSVD) ;
    
    [ZZ,YY,WW] = SVDT(AA'*BB) ;
    YY
    
    % THEREFORE, THE COMPLEMENTARY MODES ARE
    % ---------------------------------------
    % % ORTHOGONAL COMPLEMENT
    % THIS IS JUST FOR VISUALIZATION PURPOSES
    PsiSEcomp = HprojSEf_operator(PsISEfBS,Mintf, Y) ;
    DATALOCaaa.TOL = 1e-8 ;
    DATALOCaaa.Mchol =Mintfinv_chol;
    [ PsiSEcomp,S,~,~] = WSVDT( PsiSEcomp,[],DATALOCaaa) ;
    
    ncomp  = size(Y,2)-size(PsISEfBS,2) ;
    PsiSEcomp = PsiSEcomp(:,1:ncomp)  ;
    
    %
    DATALOC.LEGEND_MODES_SE = 'Complementary' ;
    %
    PlotModesSE_SubdomainLevel(DATALOC,[],PsiSEcomp,MESH) ;
    
end



disp(['------------------------------------------------'])
disp(['DETERMINATION OF BUBBLE MODES (IF ANY)'])
disp(['------------------------------------------------'])

disp(['Principal angles formed by REACTIVE FORCES and UpsilonDEF(b,:)']) ;

% PsiSEf*inv(Mintf)*PsiSEf = ident
% Let us now orthogonalize UpsilonDEF(b,:) as well, but with respect to
% Mintf
DATALOCddd = [];
[E,S_h,J ]= WSVDT( BasisUdeform.UpsilonDEF(b,:),Mintf,DATALOCddd) ;
[C,Sigma,L] = SVDT(Y'*E) ;
alpha = real(acosd(Sigma))
DATAoffline = DefaultField(DATAoffline,'MAXIMUM_ANGLE_DEF_REACT_FOR_STRAIN_MODES',89) ; % Degrees
TOL_ANGL = DATAoffline.MAXIMUM_ANGLE_DEF_REACT_FOR_STRAIN_MODES ;
disp(['Def. modes whose angle is above =',num2str(TOL_ANGL),' degrees will be considered bubble modes'])
r = find(alpha < TOL_ANGL ) ;
i_bubble = setdiff(1:length(alpha),r) ;

disp(['Number of coarse-scale strain/stresses =',num2str(length(r))]) ;
disp(['Number of bubble modes =',num2str(length(i_bubble))]) ;

% STRAIN MODES
PhiDEF = BasisUdeform.UpsilonDEF*J*L(:,r) ;

DATALOC.NAME_BASE = 'STRAIN';
PlotModesDEF_SubdomainLevel(DATALOC,PhiDEF,MESH);



% BUBBLE MODES
GammaBUB = BasisUdeform.UpsilonDEF*J*L(:,i_bubble) ;


DATALOC.NAME_BASE = 'BUBBLE';
PlotModesDEF_SubdomainLevel(DATALOC,GammaBUB,MESH);

% EFFECTIVE SELF-EQUILIBRATED MODES
PsiSEf = Y*C(:,r) ;
DATALOC.LEGEND_MODES_SE = 'EFFECTIVE_SE' ;
%
PlotModesSE_SubdomainLevel(DATALOC,[],PsiSEf,MESH) ;

% NON-EFFECTIVE SELF-EQUILIBRATED MODES
PsiSEf_non = Y*C(:,i_bubble) ;
DATALOC.LEGEND_MODES_SE = 'NONEFFECTIVE_SE' ;
%
PlotModesSE_SubdomainLevel(DATALOC,[],PsiSEf_non,MESH) ;




