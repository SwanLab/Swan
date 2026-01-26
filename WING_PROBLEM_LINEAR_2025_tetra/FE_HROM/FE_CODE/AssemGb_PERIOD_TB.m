function [G] = AssemGb_PERIOD_TB(nDOFr,nods,nods_w,nDOFm,nodm,nodm_w)
% See Plates_comp_homogenization.pdf
% Assembly of Gb matrix 
% -----------------------
DOFr = zeros(nDOFr,1) ;  % List of slave nodes
dR = zeros(nDOFr,1) ;    % Prescribed displacements
% -------------------------------------
%  Slave DOFs    associated to u    
% ---------------------------------------
iACUM = 1; n = length(nods) ;  
idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods-1)+idof];
dR(ROWS)   = uB(idof,:)'        ;
iACUM = iACUM+n-1 ;  




% ---------------------------------- 
% SLAVES DOFs associated to v    
% -------------------------------------
iACUM = iACUM + 1;  n = length(nods) ; 
idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods-1)+idof];
dR(ROWS)   = uB(idof,:)'        ;
iACUM = iACUM+n-1 ;  
% SLAVES DOFs associated to v    
% --------------------------------------
iACUM = iACUM + 1;  n = length(nods_w) ; 
idof = 3 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods_w-1)+idof];
dR(ROWS)   = uB_w(idof,:)'        ;


% 
% % ---------------------------------------
% iACUM = 1; n = length(nodm) ;  
% idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
% DOFm(ROWS) = [  3*(nodm-1)+idof];
% iACUM = iACUM+n-1 ;  
% % ---------------------------------- 
% % SLAVES DOFm associated to v    
% % -------------------------------------
% iACUM = iACUM + 1;  n = length(nodm) ; 
% idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
% DOFm(ROWS) = [  3*(nodm-1)+idof];
% iACUM = iACUM+n-1 ;  
% % SLAVES DOFm associated to v    
% % --------------------------------------
% iACUM = iACUM + 1;  n = length(nodm_w) ; 
% idof = 3 ; ROWS = iACUM:iACUM+n-1 ; 
% DOFm(ROWS) = [  3*(nodm_w-1)+idof]; 
%  
%  
%  