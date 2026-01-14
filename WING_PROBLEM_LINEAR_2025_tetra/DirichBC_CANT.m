function [DIRICHLET] = DirichBC_CANT(DATAINPUT)

% % -----------------------------------------------------------
% % 4. Dirichlet boundary conditions (prescribed displacements)
% % -----------------------------------------------------------
tEND = DATAINPUT.DATAcommon.tEND;
icond = 1; % Number of condition (fixed along time )
%---------------------------------------------------------------
%TypeFunctionDisplacementTRAINING = 'RIGID_BODY_WARPING_ON_SURFACES'; %
%ROTATION = DATAINPUT.PARAM ; 

SURFACES = [1] ; 
DIRICHLET = [] ; 

for icond = 1:length(SURFACES)

DIRICHLET(icond).NUMBER_SURFACE = SURFACES(icond);   %  
iloadstate = 1;  
DIRICHLET(icond).PRESCRIBED_DISP = [] ;
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE = {0,0,0,0,0,0,0} ;  % %  
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  %  
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [0,tEND];  % 

end
 