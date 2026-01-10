function [DIRICHLET] = DirichBound_CELL_MULT2d(t0,tEND,AMPLITUDE,SURFACES)

% % -----------------------------------------------------------
% % 4. Dirichlet boundary conditions (prescribed displacements)
% % -----------------------------------------------------------
for icond = 1:length(SURFACES)
%---------------------------------------------------------------
DIRICHLET(icond).NUMBER_SURFACE = SURFACES(icond) ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
iloadstate = 1;  % Loading stage
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE = AMPLITUDE(:,icond) ;  % % Translation (centroid face)
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [t0,tEND];

end

% 
% %%%%%%
% 
% icond = 2; % Number of condition
% %---------------------------------------------------------------
% DIRICHLET(icond).NUMBER_SURFACE = SURFACES(2) ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
% iloadstate = 1;  % Loading stage
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE =[0,0,0] ;  % % Translation (centroid face)
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [t0,tEND];  %
% %DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).ISLOCAL = 0;   
% 
% 
% icond = 3; % Number of condition
% %---------------------------------------------------------------
% DIRICHLET(icond).NUMBER_SURFACE = SURFACES(3) ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
% iloadstate = 1;  % Loading stage
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE =[0,0,0]  ;  % % Translation (centroid face)
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [t0,tEND];  %
% %DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).ISLOCAL = 0;   
% 
% 
% 
% icond = 4; % Number of condition
% %---------------------------------------------------------------
% DIRICHLET(icond).NUMBER_SURFACE = SURFACES(4) ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
% iloadstate = 1;  % Loading stage
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE =[0,0,0]  ;  % % Translation (centroid face)
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
% DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [t0,tEND];  %
% %DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).ISLOCAL = 0;   
