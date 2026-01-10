function PhiFLUC = FilterFluctuationModes(DATAcommon,PsiRBf,PhiFLUC,PsiDEFf,Vrb,Mintf)


    % WORK DONE BY THE RESULTANT FORCES
    WorkResForces = PsiRBf'*PhiFLUC ;
    disp(['NORM WORK DONE BY  RESULTANT FORCES over fluctuation displacements =',num2str(norm(WorkResForces,'fro'))]) ;
    % WORK DONE BY THE self-equilibrated  FORCES
    % COMPUTING PRINCIPAL ANGLES BETWEEN RESULTANT FORCES
    DATALOC.RELATIVE_SVD = 1 ;
    TOL_SVD  = 1e-6 ;
    [QA,SA,~] = SVDT(PsiDEFf,TOL_SVD,DATALOC) ;
    [QB_fl,SB,~] = SVDT(PhiFLUC,TOL_SVD,DATALOC) ;
    [Yfluc,ScosineFLUC,Z] = SVDT(QB_fl'*QA) ;
    PrincipalAnglesFluct = acosd(ScosineFLUC) ;
 
    [QB_rb,SB,~] = SVDT(Vrb,TOL_SVD,DATALOC) ;
    [Yrb,ScosineRB,Z] = SVDT(QA'*QB_rb) ;
    PrincipalAnglesRB = acosd(ScosineRB) ;
     disp('Principal angles fluctuations')
    PrincipalAnglesRB
     
        disp('Principal angles fluctuations')
    PrincipalAnglesFluct
    
     
    switch DATAcommon.RESTRICT_NUMBER_FLUCTUATION_MODES_CRITERION
        case 'BY_ANGLES_Vrb'
            
            LLL = find(PrincipalAnglesFluct<= PrincipalAnglesRB(end)) ;
             PhiFLUC = QB_fl*Yfluc(:,LLL) ;
             
             [PhiFLUC,~,~] = WSVDT(PhiFLUC,Mintf)  ; 
             
        otherwise
            
            
    end