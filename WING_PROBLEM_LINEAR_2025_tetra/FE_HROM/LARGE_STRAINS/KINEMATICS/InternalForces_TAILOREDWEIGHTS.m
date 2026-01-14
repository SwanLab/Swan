function  Fint = InternalForces_TAILOREDWEIGHTS(OPERFE,PoneST,PK2STRESS,DATA) ;  
% Internal forces, large strains, in terms of the Bst, W operators, as well
% as the 1st PK stress (stacked) 
% Tailored weights for the ECM, see 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/106_ECMtailoredTOmodes/01_VIABILITY.mlx
% 
if nargin == 0
    load('tmp.mat')
end


if ~isempty(PoneST)
    
    
    Fint = OPERFE.BstW'*PoneST ;
else
   
    
    Fint = OPERFE.BstW'*PK2STRESS ;
end

