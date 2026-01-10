function  [PhiDEF,PsiDEFf,MdomCHOL] =   DeformationalModes_EIFE(faceDOFSall,AsnapDISP,PsiDEFf,PhiRB,Mdom)  ;
% Computation of deformational modes, EIFE METHOD 
% Adaptation of DeformationalModes_1dom, developed in 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
%  
% New scenario explained in 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx
% 
if nargin == 0
    load('tmp.mat')
    
end

% STEP 1 ) DEFORMATIONAL COMPONENT OF THE DISPLACEMENT SNAPSHOTS 
AsnapDISPdef =  SprojDEF_operator(PhiRB,Mdom,AsnapDISP) ;
f = faceDOFSall ;
% STEP 2
% SINGULAR VALUE DECOMPOSITION OF THE RESTRICTION OF THE DEFORMATIONAL
% DISPLACEMENTS TO THE BOUNDARY OF THE EIF ELEMENT  
[U_a,S_a,V_a ]= SVDT(AsnapDISPdef(f,:),0) ;  
% Coefficients 
c_a = bsxfun(@times,V_a',1./S_a)' ;   % These are the coefficients of the projection onto U_a
% In principle, one could make    PhiDEF  =  AsnapDISPdef*c ; 
% YEt in proceeding this way,  we do not ensure one of the conditions upon
% which the proposed approach is founded, namely, that PsiDEFf'*PhiDEFf
% should be square (and invertible)
% To ensure this condition, we must seek the subspace of the column space
% of U_a that is "effective", that is, that does work; this is given by  
% U_a*COEFF, where 
COEFF = U_a'*PsiDEFf ;
[UU,SS,VV] = SVDT(COEFF);  % SVD coefficients
TOL_illposedness = 1e-4 ;  % This may be raised if required
if find(SS/SS(1) <TOL_illposedness)
    error('MAtrix Hdef = PsiDEFf^T*PhiDEFf might result rank deficient... REdesign the training process')
end
[PhiDEF,~,~,MdomCHOL] = WSVDT(AsnapDISPdef*(c_a*COEFF),Mdom) ;

 