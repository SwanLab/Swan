function [DOFm ROWSmst] = AssemDOFmst_PERIOD_TB(nDOFm,nodm,nodm_w)
% See Plates_comp_homogenization.pdf
% Assembly of master DOFs 
% -----------------------
DOFm = zeros(nDOFm,1) ;  % List of slave nodes
% -------------------------------------
%  Slave DOFm    associated to u    
% ---------------------------------------
iACUM = 1; n = length(nodm) ;  
idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
DOFm(ROWS) = [  3*(nodm-1)+idof];
ROWSmst.u = ROWS ; 
iACUM = iACUM+n-1 ;  
% ---------------------------------- 
% SLAVES DOFm associated to v    
% -------------------------------------
iACUM = iACUM + 1;  n = length(nodm) ; 
idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
DOFm(ROWS) = [  3*(nodm-1)+idof];
ROWSmst.v = ROWS ; 
iACUM = iACUM+n-1 ;  
% SLAVES DOFm associated to v    
% --------------------------------------
iACUM = iACUM + 1;  n = length(nodm_w) ; 
idof = 3 ; ROWS = iACUM:iACUM+n-1 ; 
ROWSmst.w = ROWS ; 
DOFm(ROWS) = [  3*(nodm_w-1)+idof]; 
 