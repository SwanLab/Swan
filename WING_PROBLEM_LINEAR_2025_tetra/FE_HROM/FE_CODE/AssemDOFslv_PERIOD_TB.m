function [DOFr dR ROWSslv] = AssemDOFslv_PERIOD_TB(nDOFr,nods,nods_w,uB,uB_w)
% See Plates_comp_homogenization.pdf
% Assembly of slave DOFs 
% -----------------------
DOFr = zeros(nDOFr,1) ;  % List of slave nodes
dR = zeros(nDOFr,1) ;    % Prescribed displacements
ROWSslv = [] ; 
% -------------------------------------
%  Slave DOFs    associated to u    
% ---------------------------------------
iACUM = 1; n = length(nods) ;  
idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods-1)+idof];
dR(ROWS)   = uB(idof,:)'        ;
ROWSslv.u = ROWS ;  
iACUM = iACUM+n-1 ;  
% ---------------------------------- 
% SLAVES DOFs associated to v    
% -------------------------------------
iACUM = iACUM + 1;  n = length(nods) ; 
idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods-1)+idof];
dR(ROWS)   = uB(idof,:)'        ;
ROWSslv.v = ROWS ;  
iACUM = iACUM+n-1 ;  
% SLAVES DOFs associated to v    
% --------------------------------------
iACUM = iACUM + 1;  n = length(nods_w) ; 
idof = 3 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods_w-1)+idof];
dR(ROWS)   = uB_w(idof,:)'        ;
 ROWSslv.w = ROWS ;  
 