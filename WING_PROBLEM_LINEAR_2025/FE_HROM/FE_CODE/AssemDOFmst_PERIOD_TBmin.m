function [DOFm ROWSmst] = AssemDOFmst_PERIOD_TBmin(nDOFm,nodm,MASTER_uvFL)
% See Plates_comp_homogenization.pdf
% Assembly of master DOFs 
% -----------------------
DOFm = zeros(nDOFm,1) ;  % List of slave nodes
% -------------------------------------
%  Slave DOFm    associated to u     (PERIODIC)
% ---------------------------------------
iACUM = 1; n = length(nodm) ;  
idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
DOFm(ROWS) = [  3*(nodm-1)+idof];
ROWSmst.u = ROWS ; 
iACUM = iACUM+n-1 ;  
% ---------------------------------- 
% SLAVES DOFm associated to v    (PERIODIC)
% -------------------------------------
iACUM = iACUM + 1;  n = length(nodm) ; 
idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
DOFm(ROWS) = [  3*(nodm-1)+idof];
ROWSmst.v = ROWS ; 
iACUM = iACUM+n-1 ;  
 % ---------------------------------- 
% SLAVES DOFm associated to w    (PERIODIC)
% -------------------------------------
iACUM = iACUM + 1;  n = length(nodm) ; 
idof = 3 ; ROWS = iACUM:iACUM+n-1 ; 
DOFm(ROWS) = [  3*(nodm-1)+idof];
ROWSmst.w = ROWS ; 
iACUM = iACUM+n-1 ;  
 % ---------------------------------- 
% SLAVES DOFm associated to u    (MINIMAL)
% -------------------------------------
iACUM = iACUM + 1;  n = length(MASTER_uvFL) ; 
idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
DOFm(ROWS) = [  3*(MASTER_uvFL-1)+idof];
ROWSmst.uMIN = ROWS ; 
iACUM = iACUM+n-1 ;  
 % ---------------------------------- 
% SLAVES DOFm associated to v    (MINIMAL)
% -------------------------------------
iACUM = iACUM + 1;  n = length(MASTER_uvFL) ; 
idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
DOFm(ROWS) = [  3*(MASTER_uvFL-1)+idof];
ROWSmst.vMIN = ROWS ; 
iACUM = iACUM+n-1 ;  
 