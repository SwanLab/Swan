function [DOFr dR ROWSslv] = AssemDOFslv_PERIOD_TBmin(nDOFr,nods,uB,SLAVE_uvFL,uBmacro)
% See Plates_comp_homogenization.pdf
% Assembly of slave DOFs 
% -----------------------
DOFr = zeros(nDOFr,1) ;  % List of slave nodes
dR = zeros(nDOFr,1) ;    % Prescribed displacements
ROWSslv = [] ; 
% -------------------------------------
%  Slave DOFs    associated to u    (PERIODIC)
% ---------------------------------------
iACUM = 1; n = length(nods) ;  
idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods-1)+idof];
dR(ROWS)   = uB(idof,:)'        ;
ROWSslv.u = ROWS ;  
iACUM = iACUM+n-1 ;  
% ---------------------------------- 
% SLAVES DOFs associated to v    (PERIODIC)
% -------------------------------------
iACUM = iACUM + 1;  n = length(nods) ; 
idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods-1)+idof];
dR(ROWS)   = uB(idof,:)'        ;
ROWSslv.v = ROWS ;  
iACUM = iACUM+n-1 ;  
% SLAVES DOFs associated to w (PERIODIC)
% --------------------------------------
iACUM = iACUM + 1;  n = length(nods) ; 
idof = 3 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(nods-1)+idof];
dR(ROWS)   = uB(idof,:)'        ;
ROWSslv.w = ROWS ;  
iACUM = iACUM+n-1 ;  
% SLAVES DOFs associated to u (minimal on top/bottom faces)
% --------------------------------------
iACUM = iACUM + 1;  n = length(SLAVE_uvFL) ; 
idof = 1 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(SLAVE_uvFL-1)+idof];
dR(ROWS)   = uBmacro.u ; 
ROWSslv.uMIN = ROWS ;  
iACUM = iACUM+n-1 ;  
% SLAVES DOFs associated to v (minimal on top/bottom faces)
% --------------------------------------
iACUM = iACUM + 1;  n = length(SLAVE_uvFL) ; 
idof = 2 ; ROWS = iACUM:iACUM+n-1 ; 
DOFr(ROWS) = [  3*(SLAVE_uvFL-1)+idof];
dR(ROWS)   = uBmacro.v ; 
ROWSslv.vMIN = ROWS ;  
iACUM = iACUM+n-1 ;  
 
 