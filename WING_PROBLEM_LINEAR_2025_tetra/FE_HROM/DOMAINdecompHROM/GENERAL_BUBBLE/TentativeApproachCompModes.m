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