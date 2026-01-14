function K = GetStiffnessMatrix_COROT_LRss(OPERFE,DATA,VAR,FgradST,celastST,LboolCallQ,D_QrotALL,KcGEOunassemGLOloc)
% Adaptation of GetStiffnessMatrix.m to the co-rotational method
% See  Small strains/Large rotations (or Small rotations/Large strains)
% JAHO, 14-feb-2025, FRIDAY, 7:43, Balmes 185,  Barcelona.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% and
%/home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf

if DATA.SMALL_STRAIN_KINEMATICS == 0
     K = KstiffLargeStrains_COROT_LRss(OPERFE,VAR.PK2STRESS,FgradST,DATA.MESH.ndim,celastST,DATA,LboolCallQ,D_QrotALL,KcGEOunassemGLOloc) ;
    
    
    
else
    error('Option not implemented')
 %   K = KstiffSmallStrains(OPERFE,FgradST,DATA.MESH.ndim,celastST,DATA) ;
end
