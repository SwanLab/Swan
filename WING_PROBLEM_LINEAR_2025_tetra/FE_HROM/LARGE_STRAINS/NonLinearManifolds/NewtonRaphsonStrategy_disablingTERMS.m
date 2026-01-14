function [VAR,CONVERGED ]= NewtonRaphsonStrategy_disablingTERMS(DATA,OPERFE,VAR,MATPRO,DOFl)

USE_GEOMETRIC_term_K = 1;  USE_WEIGHT_TERM_k = 1;
[VAR_new,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMDw(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,USE_WEIGHT_TERM_k) ;

if CONVERGED == 0 && DATA.DisableGeometricAndWeightTermWhenNoConvergence
    [VAR,CONVERGED] = Hier_NewtonRapshonStaticLarge_manifoldMDw(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
        USE_WEIGHT_TERM_k) ;
    
    %             elseif CONVERGED == 0 && DATA.USE_BACKTRACK
    %                 % https://chatgpt.com/share/69013953-5be8-8013-ac25-618113b6dbf2
    %                 error('Option not operative')
    %                  [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMDw_CGPT(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
    %                 USE_WEIGHT_TERM_k) ;
else
    
    VAR = VAR_new;
end