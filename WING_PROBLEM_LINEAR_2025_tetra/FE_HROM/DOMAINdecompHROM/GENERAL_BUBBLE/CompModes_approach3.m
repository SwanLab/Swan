    disp(['Method 2 decomposition "inelastic" modes. Principal angles PsiSE and UpsilonDEF --> GammaBUB'])
    % See  %
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/05_Assessment.mlx
    % -------------------------------------------------------------------------------------------------------------
    % UpsilonDEF = [PhiDEF,BasisUdeform.PhiDEFcomp] ;
    % [Z,SS,VV] =SVDT(UpsilonDEF(b,:)) ;
    %
    % Y =   SVDT([PsiSEf,PhiDEF(b,:)]) ;
    
    UpsilonDEF = [BasisUdeform.PhiDEFcomp] ;
    [Z,SS,VV] =SVDT(UpsilonDEF(b,:)) ;
    
    Y =   SVDT([PsiSEf]) ;
    
    
    
    disp('-------------------------------------------------')
    disp('PRINCIPAL  ANGLES between PsiSEf    AND   Z)')
    disp('-------------------------------------------------')
    % [ZZ,Sz,Vz]  = SVDT(Z) ;
    % [YY,Sy,Vy]  = SVDT(Y) ;
    
    
    [uu,ss,vv] = SVDT(Z'*Y) ;
    theta = real(acosd(ss))
    
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
    
    %
    