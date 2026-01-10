
    disp('----------------------------------------------------------------------')
    disp(['Computing self-equilibrated modes counterpart of the deformational modes '])
    disp([' via right-singular vectors of displacements'])
    disp(['Just for completeness, because the "inelastic" SE modes are not used in this approach '])
    
    disp('----------------------------------------------------------------------')
    % A = U*S*V^T -->   U = (A*V)*inv(S)
    
    disp(['This operation may be optimized...'])
    PsiSEf_viaCOEFF = cell2mat(AsnapREACse)*BasisUdeform.V_upsilonDEF ;
    % Normalization
    DATALOCcc.Mchol = [] ;
    [ Y,~,~,~] = WSVDT( PsiSEf_viaCOEFF,MintfINV,DATALOCcc) ;
    
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
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES    BY PsiSEcomp AND  PhiDEFb  (= BasisUdeform.PhiDEFbs(b,:))')
    disp('-------------------------------------------------')
    
    [PhiDEFb,SSS,VVV] =WSVDT(BasisUdeform.PhiDEFbs(b,:),Mintf,[]) ;
    [uu,ss,vv] = SVDT(PsiSEcomp'*PhiDEFb) ;
    alpha = real(acosd(ss))