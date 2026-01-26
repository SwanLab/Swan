function  [PsiSEf,PhiDEF,GammaBUB] =   SelfEquilibratedModes_bub(PsiRBf,Mintf,AsnapREAC,DATALOC,DATAoffline,MESH,BasisUdeform,Mdom)
% COMPUTATION OF SELF-EQUILIBRATED MODES
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/05_Assessment.mlx
% JAHO, 17-Oct-2023

if nargin == 0
    load('tmp1.mat')
end
disp(['-----------------------------------------------------------'])
disp(['Determination of self-equilibrated modes'])
disp(['-----------------------------------------------------------'])

AsnapREACse = cell(size(AsnapREAC)) ;
b = MESH.faceDOFSall ;  % Boundary DOFs
% Self-equilibrated part (in principle, not necessary if one trains without external forces, just prescribed displacements)
for iproj = 1:length(AsnapREAC)
    [AsnapREACse{iproj},MintfINV] =  HprojSEf_operator(PsiRBf,Mintf,AsnapREAC{iproj}(b,:)) ;
end
SNAPbasic = cell2mat(AsnapREACse(DATAoffline.INDEX_BASIC_TESTS)) ;  % Elastic
disp('Computing  "basic" self-equilibrated modes')
DATALOCww.TOL = 1e-8;
%DATALOCww.Mchol = Mintfinv_chol ;
[ PsISEfBS,Sbs,~,Mintfinv_chol] = WSVDT( SNAPbasic,MintfINV,DATALOCww) ;

DATALOC.LEGEND_MODES_SE = 'Basic' ;
PlotModesSE_SubdomainLevel(DATALOC,PsiRBf,PsISEfBS,MESH) ;

disp(['Approach in which Self-equilibrated modes = Basis Modes'])

PsiSEf = PsISEfBS ;
% Like wise
PhiDEF = BasisUdeform.PhiDEFbs ;
% ------------------------------------------------------------------


CHECK_PRINCIPAL_ANGLES_ORTHOGONAL_COMPL_MODES = 0;
if CHECK_PRINCIPAL_ANGLES_ORTHOGONAL_COMPL_MODES == 1
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES    BY PsiSEf = PsISEfBS    AND   BasisUdeform.PhiDEFcomp(b,:)  (= BasisUdeform.PhiDEFbs(b,:))')
    disp('-------------------------------------------------')
    [PhiDEFcompB,SSS,VVV] =WSVDT(BasisUdeform.PhiDEFcomp(b,:),Mintf,[]) ;
    [uu,ss,vv] = SVDT(PsiSEf'*PhiDEFcompB) ;
    beta = real(acosd(ss))
end
CHECK_COMPLEMENTARY_MODES = 0 ; % Just if one wants to see how the discarded reactive modes look like
if CHECK_COMPLEMENTARY_MODES == 1
    ComplModesMetho1; % See inside. One of the approaches tested before arriving at the "good" one
end


%
METHOD_nullspace = 0 ;

if   METHOD_nullspace == 1
    disp(['Method based on the kernel of PsiSE^T'])
    [Kker,SS,VV] = svd(PsiSEf) ;
    Kker = Kker(:,size(PsiSEf,2)+1:end) ;
    UpsilonDEF = [PhiDEF,BasisUdeform.PhiDEFcomp] ;
    [Z,SS,VV] =SVDT(UpsilonDEF(b,:)) ;  % Z*SS*VV^T = Ub  --> U(b,:)*VV*inv(SS) = Z
    vvINCS = bsxfun(@times,VV',1./SS)' ;
    [uu,ss,vv] = SVDT(Z'*Kker) ;
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES between kernel(PsiSE^T)    AND   UpsilonDEF(b,:))')
    disp('-------------------------------------------------')
    theta = real(acosd(ss))
    ninel = size(BasisUdeform.PhiDEFcomp,2) ;
    GammaBUB = UpsilonDEF*vvINCS*uu(:,1:ninel) ;
    [GammaBUB,SS,VV] = WSVDT(GammaBUB,Mdom) ;
    disp('Checking that GammaBUB and PhiDEFf are linearly independent')
    disp('If this is so, all the angles formed by the corresponding subspaces should greater than 0') ;
    [uuu,sss,vvv] = SVDT(GammaBUB'*Mdom*PhiDEF) ;
    alpha = real(acosd(sss))
    % Let
    if any(alpha < 1.0)  %%%The
        warning('Some bubble modes are almost linearly dependent to the basic modes')
    end
    disp(['Ckecking that the modes are indeed free-residual bubble modes '])
    [Z,SS,VV] =WSVDT(GammaBUB(b,:),Mintf) ;  % Z*SS*VV^T = Ub  --> U(b,:)*VV*inv(SS) = Z
    [uu,ss,vv] = SVDT(Z'*PsiSEf) ;
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES between  PsiSE     AND  GammaBUB(b,:) )')
    disp('-------------------------------------------------')
    theta = real(acosd(ss))
    DATALOC.NAME_BASE = 'BUBBLE';
    PlotModesDEF_SubdomainLevel(DATALOC,GammaBUB,MESH);
end



METHOD_nullspace2 = 1 ;

if   METHOD_nullspace2 == 1
    disp(['Method based on orthogonal projections  '])
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/06_ForcingBubbles.mlx
    %
    UpsilonDEF = [PhiDEF,BasisUdeform.PhiDEFcomp] ;
    epsilon_d = 1e-6;
    DATALOCaaa.RELATIVE_SVD = 1;
    [Dbar,S_d,V_d] = SVDT(UpsilonDEF(b,:),epsilon_d,DATALOC) ;
    VSd = bsxfun(@times,V_d',1./S_d)' ;
    D = UpsilonDEF*VSd ;
    
    if size(D,2) ~= size(UpsilonDEF,2)
        error('To be implemented')
    end
    p  = size(D,2) - size(PhiDEF,2) ;
    [P,SS1,VV] = SVDT(PsiSEf) ;
    % Projectio of Dbar onto the orthogonal complement of span(P)
    T = Dbar- P*(P'*Dbar);
    [G,St,Vt] = SVDT(T) ;
    [Jbar, SJ,VJ] = SVDT(G'*Dbar) ;
    disp(['Principal Angles formed by the kernel of PsiSE and the columnspace of UpsilonDEF'])
    disp(['If the first ',num2str(p),' are zero, it means the bubble condition is exactly met (no need for approximations)'])
    real(acosd(SJ))
    Y = D*VJ(:,1:p) ;
    
    s = 1:size(Y,1) ;
    s(b)  = [] ;
    GammaBUB = zeros(size(Y,1),p) ;
    GammaBUB(s,:) = Y(s,:) ;
    GammaBUB(b,:) = G*Jbar(:,1:p);
    [GammaBUB,SS,VV] = WSVDT(GammaBUB,Mdom) ;
    
    
    disp('Checking that GammaBUB and PhiDEFf are linearly independent')
    disp('If this is so, all the angles formed by the corresponding subspaces should greater than 0') ;
    [uuu,sss,vvv] = SVDT(GammaBUB'*Mdom*PhiDEF) ;
    alpha = real(acosd(sss))
    % Let
    if any(alpha < 1.0)  %%%The
        warning('Some bubble modes are almost linearly dependent to the basic modes')
    end
    disp(['Ckecking that the modes are indeed free-residual bubble modes '])
    [Z,SS,VV] =WSVDT(GammaBUB(b,:),Mintf) ;  % Z*SS*VV^T = Ub  --> U(b,:)*VV*inv(SS) = Z
    [uu,ss,vv] = SVDT(Z'*PsiSEf) ;
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES between  PsiSE     AND  GammaBUB(b,:) )')
    disp('-------------------------------------------------')
    theta = real(acosd(ss))
    DATALOC.NAME_BASE = 'BUBBLE';
    PlotModesDEF_SubdomainLevel(DATALOC,GammaBUB,MESH);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%% TESTED APPROACHES  (NOT SATISFACTORY IN THE END )%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
METHOD_2 = 0 ;
if METHOD_2 == 1
    CompModes_approach3.m; % Another tentative approach, see inside
end

TENTATIVE_APPROACH_DIDNOTWORK = 0;

if TENTATIVE_APPROACH_DIDNOTWORK ==1
    TentativeApproachCompModes;
end

METHOD_1 = 0 ;
if METHOD_1 == 1
    ComlModesOblique; % This yet another approaches we tested before arriving at the final gone, see below
    
end



