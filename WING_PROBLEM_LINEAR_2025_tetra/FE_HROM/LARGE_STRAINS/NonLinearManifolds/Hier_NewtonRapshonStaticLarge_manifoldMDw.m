function [VAR,CONVERGED] = Hier_NewtonRapshonStaticLarge_manifoldMDw(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
    USE_WEIGHT_TERM_k)

disp('Disabling geometric matrix')
USE_GEOMETRIC_term_K = 0 ;   USE_WEIGHT_TERM_k = 1;
[VAR_new,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMDw(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
    USE_WEIGHT_TERM_k) ;

if  CONVERGED == 0
    disp('Disabling  weight matrix (but not geometric one)')
    USE_GEOMETRIC_term_K = 1 ;   USE_WEIGHT_TERM_k = 0;
    [VAR_new,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMDw(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
        USE_WEIGHT_TERM_k) ;
    
    if  CONVERGED == 0
        disp('Disabling both geometric and weight matrix  ')
        USE_GEOMETRIC_term_K = 0 ;   USE_WEIGHT_TERM_k = 0;
        [VAR,CONVERGED ]= NewtonRapshonStaticLarge_manifoldMDw(DATA,OPERFE,VAR,MATPRO,DOFl,USE_GEOMETRIC_term_K,...
            USE_WEIGHT_TERM_k) ;
    
    else
        VAR = VAR_new;
    end
else
    VAR = VAR_new;
end
