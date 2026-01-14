function [StwoST,celasST,detFgrad]= NeoHookStressStrain(detFgrad,EgreenlST,MATPRO,DATA)

if nargin == 0
    load('tmp.mat')
    %
    F = eye(3) ;
    F(1,1) = [1.0001] ;
    F(2,1) = 0.0001 ;
    detF = det(F) ;
    E = 0.5*(F'*F-eye(3)) ;
    EgreenlST = [E(1,1)
        E(2,2)
        2*E(1,2) ] ;
    
    
    %             F_mat = [F1  F6 F5
    %          F9  F2 F4
    %         F8   F7 F3 ] ;
    FgradST = [F(1,1)
        F(2,2)
        F(1,2)
        F(2,1)
        
        ] ;
    DATA.MESH.nstrain  =3;
    DATA.MESH.ndim  =2 ;
    
    C =F'*F ;
    inv(C)
    
    Celas =          [91731.2661498708          37467.7002583979                         0
        37467.7002583979          91731.2661498708                         0
        0                         0          27131.7829457364] ;
    
    MATPRO.mu_0 = MATPRO.mu_0(1) ;
    MATPRO.lambda_0 = MATPRO.lambda_0(1) ;
end


% Inverse of F^T*F, given as data the determinant of (F^T*F), which is
% detFgrad.^2
CrightINV = CrightINV_fromEgreenlST(EgreenlST,DATA.MESH.nstrain,detFgrad.^2) ;
%  PK2 stresses as a function  CrightINV
% \voigt{\Stwo}_a = \lambda_0 \logn( \deter{\,\Fgrad} )     {\CrightINV}_a  + \mu_0 (\IDENTITY{a} -  {\CrightINV}_a)

% Let's start with the term     (\IDENTITY{A} -  {\CrightINV}_a)
StwoST  = -CrightINV ;
for idim = 1:DATA.MESH.ndim
    LOCROWS = idim:DATA.MESH.nstrain:length(StwoST) ;
    StwoST(LOCROWS) = 1+StwoST(LOCROWS) ;
end

% \lambda_0 \logn( \deter{\,\Fgrad} )     {\CrightINV}_a  + \mu_0  (******)

for istrain = 1:DATA.MESH.nstrain
    ROWS = istrain:DATA.MESH.nstrain:length(StwoST) ;
    StwoST(ROWS) = MATPRO.mu_0.*StwoST(ROWS) + MATPRO.lambda_0.*log(detFgrad).*CrightINV(ROWS) ;
end


%%%% Tangent matrix
% ***********************
% \begin{equation}
%  \ccelas{A}{B}{C}{D}  = \lambda_0 \overbrace{\CCrightINV{A}{B} \CCrightINV{C}{D}}^{\ddSS \IIvol{A}{B}{C}{D}}
% +  2\mu \overbrace{\Par{ \dfrac{1}{2}  (\CCrightINV{A}{C} \CCrightINV{B}{D}  + \CCrightINV{A}{D} \CCrightINV{B}{C}) }}^{\ddSS \IIsym{A}{B}{C}{D}}
% \end{equation}
%
% where
%
% \begin{equation}
%  \mu \defeq \mu_0 -  \lambda_0 \logn( \deter{\,\Fgrad} )
% \end{equation}

mu = MATPRO.mu_0 -MATPRO.lambda_0.*log(detFgrad) ;
% Volumetric component  % \CCrightINV{A}{B} \CCrightINV{C}{D}
celasST = Ivol_LargeStrains(CrightINV,DATA.MESH.nstrain) ;
% Identity component  %  \dfrac{1}{2}  (\CCrightINV{A}{C} \CCrightINV{B}{D}  + \CCrightINV{A}{D} \CCrightINV{B}{C})
Isym = Isym_LargeStrains(CrightINV,DATA.MESH.nstrain) ;

for istrain = 1:DATA.MESH.nstrain
    ROWS = istrain:DATA.MESH.nstrain:length(StwoST) ;
    celasST(ROWS,:) =  MATPRO.lambda_0.*celasST(ROWS,:)    + 2*mu.*Isym(ROWS,:)  ;
end

