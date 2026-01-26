function K = GetStiffnessMatrix_bub(OPERFE,DATA,VAR,FgradST,celastST)
 
    if DATA.SMALL_STRAIN_KINEMATICS == 0
        error('Option not implemented yet ')
        K = KstiffLargeStrains(OPERFE,VAR.PK2STRESS,FgradST,DATA.MESH.ndim,celastST) ;
    else
        K = KstiffSmallStrains_bub(OPERFE,FgradST,DATA.MESH.ndim,celastST) ;
    end
 