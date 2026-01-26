disp(['Method 1 decomposition "inelastic" modes. Oblique projection rather than orthogonal projection --> GammaBUB'])
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/05_Assessment.mlx
    cOB = (PhiDEF(b,:)'*Mintf*PhiDEF(b,:))\(PhiDEF(b,:)'*Mintf*BasisUdeform.PhiDEFcomp(b,:)) ;
    PhiDEFcomp = BasisUdeform.PhiDEFcomp - PhiDEF*cOB ;
    [GammaBUB,SSS,VVV] =WSVDT(PhiDEFcomp,Mdom,[]) ;
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES    BY PsiSEf = PsISEfBS    AND   GammaBUB)')
    disp('-------------------------------------------------')
    [GammaBUBb,SSS,VVV] =WSVDT(GammaBUB(b,:),Mintf,[]) ;
    [uu,ss,vv] = SVDT(GammaBUBb'*PsiSEf) ;
    gamma = real(acosd(ss))
    
    
    % PRINCIPAL ANGLES COMPUTED ALTERNATIVELY BY USING THE STANDARD DEFINITION
    [GammaBUBb_s,SSS2,VV2V] =SVDT(GammaBUBb) ;
    [PsiSEf_s,SSS3,VVV3] =SVDT(PsiSEf) ;
    [uu,ss,vv] = SVDT(GammaBUBb_s'*PsiSEf_s) ;
    gamma = real(acosd(ss))
    
    
    disp(['Deformed shaped bubble modes ---ordered according to the angles formed with PsiSEf'])
    GammaBUBsorted = GammaBUB*uu ;
    
    DATALOC.NAME_BASE = 'BUBBLE_SORTED';
    PlotModesDEF_SubdomainLevel(DATALOC,GammaBUBsorted,MESH);
    
    for imodeINEL = 1:size(GammaBUBsorted,2)
        disp(['Quasi-bubble mode GammaBUBsorted = ',num2str(imodeINEL)])
        normB = sqrt(GammaBUBsorted(b,imodeINEL)'*Mdom(b,b)*GammaBUBsorted(b,imodeINEL)) ;
        normT =     sqrt(GammaBUBsorted(:,imodeINEL)'*Mdom*GammaBUBsorted(:,imodeINEL)) ;
        %  % This is one
        disp(['Norm boundary displacements   (over 1)',num2str(normB/normT)])
    end