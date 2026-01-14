function  [PsiSEf,PhiDEF,GammaBUB] =   BubbleModes_BasicSE(PsiRBf,Mintf,AsnapREAC,DATALOC,DATAoffline,MESH,BasisUdeform,Mdom)
% COMPUTATION OF SELF-EQUILIBRATED MODES
%/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/06_ForcingBubbles.mlx
% JAHO, 23-Oct-2023
% -------------------------------------------
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

%
INDEX_FUNDAMENTAL_MODES  =  DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.INDEX_FUNDAMENTAL_MODES ;

disp('Computing  "basic" self-equilibrated modes')
DATALOCww.TOL = 1e-8;
if isempty(INDEX_FUNDAMENTAL_MODES)
    
    
    %DATALOCww.Mchol = Mintfinv_chol ;
    [ PsISEfBS,Sbs,~,Mintfinv_chol] = WSVDT( SNAPbasic,MintfINV,DATALOCww) ;
    
else
    %  See
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/01_BEAMQ8_ELAST.mlx
    % Separated treatment of "fundamental" basic modes (for instance, Q8 element for capturing beam bending behavior)
    % See also /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/DeformModesNONL.m
    [ PsISEfBS_fund,Sbs_fund,~,Mintfinv_chol] = WSVDT( SNAPbasic(:,INDEX_FUNDAMENTAL_MODES),MintfINV,DATALOCww) ;
    INDEX_remaining_modes = setdiff(1:size(SNAPbasic,2),INDEX_FUNDAMENTAL_MODES) ; 
    DATALOCww.Mchol = Mintfinv_chol ; 
     [ PsISEfBS_remaining,Sbsremain,~,Mintfinv_chol] = WSVDT( SNAPbasic(:,INDEX_remaining_modes),[],DATALOCww) ;
    
    PsISEfBS = [PsISEfBS_fund,PsISEfBS_remaining] ; 
    
     
    
    if size(PsISEfBS,2) ~= DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_BASIC_DEFORMATIONAL_MODES
        error(['The number of   elastic modes should be equal to ',num2str(DATAoffline.INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_BASIC_DEFORMATIONAL_MODES)])
    end
end






DATALOC.LEGEND_MODES_SE = 'Basic' ;
PlotModesSE_SubdomainLevel(DATALOC,PsiRBf,PsISEfBS,MESH) ;

disp(['RECALL THAT THE APPROACH USED IS THAT IN WHICH THE  ONLY Self-equilibrated modes IN THE FORMULATION ARE BasiC Modes (ELASTIC)'])

PsiSEf = PsISEfBS ;
% Like wise
PhiDEF = BasisUdeform.PhiDEFbs ;
% ------------------------------------------------------------------

CHECK_MATRIX_H_BASIC_MODES  = 1 ; 
disp('------------------------------------------------------------------------------------------------')
if  CHECK_MATRIX_H_BASIC_MODES == 1
    % Principal angles formed by PhiDEFb and PsiSEf
    PhiDEFb = SVDT(PhiDEF(b,:)) ; 
    PsiSEf_orth = SVDT(PsiSEf) ; 
    [UUU,PrincipalAngles_cosines,VVV] = SVDT(PsiSEf_orth'*PhiDEFb) ; 
    disp(['Cosine principal angles  formed by basic deformational and self-equilibrated modes'])
      PrincipalAngles_cosines 
      TOL_limit_cosines = 1e-4 ; 
     if PrincipalAngles_cosines(end) < TOL_limit_cosines
         warning(['Smallest cosine principal angle is below '],num2str(TOL_limit_cosines))
         disp(['This might be conducive to ill-posed problems'])
         pause
     
     end
    
    
    disp('------------------------------------------------------------------------------------------------')

    
end


CHECK_PRINCIPAL_ANGLES_ORTHOGONAL_COMPL_MODES = 0;
if CHECK_PRINCIPAL_ANGLES_ORTHOGONAL_COMPL_MODES == 1
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES    BY PsiSEf = PsISEfBS    AND   BasisUdeform.PhiDEFcomp(b,:)  (= BasisUdeform.PhiDEFbs(b,:))')
    disp('-------------------------------------------------')
    [PhiDEFcompB,SSS,VVV] =WSVDT(BasisUdeform.PhiDEFcomp(b,:),Mintf,[]) ;
    [uu,ss,vv] = SVDT(PsiSEf'*PhiDEFcompB) ;
    beta = real(acosd(ss))
end

%


disp(['Method based on orthogonal projections  '])
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/06_ForcingBubbles.mlx
%
UpsilonDEF = [PhiDEF,BasisUdeform.PhiDEFcomp] ;
epsilon_d = 1e-6;
DATALOC.RELATIVE_SVD = 1;
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
