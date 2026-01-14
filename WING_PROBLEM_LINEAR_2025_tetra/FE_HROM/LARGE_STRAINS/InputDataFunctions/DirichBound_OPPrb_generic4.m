function [DIRICHLET] = DirichBound_OPPrb_generic4(t0,tEND,AMPLITUDE,surf_to_move,rest_surf)

% % -----------------------------------------------------------
% % 4. Dirichlet boundary conditions (prescribed displacements)
% % -----------------------------------------------------------
icond = 1; % Number of condition
%---------------------------------------------------------------
DIRICHLET(icond).NUMBER_SURFACE = surf_to_move ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
iloadstate = 1;  % Loading stage
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE =AMPLITUDE(:)' ;  % % Translation (centroid face)
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [t0,tEND];

%%%%%%

for icondLOC = 1:length(rest_surf)
icond = 1+icondLOC; % Number of condition
%---------------------------------------------------------------
DIRICHLET(icond).NUMBER_SURFACE = rest_surf(icondLOC) ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
iloadstate = 1;  % Loading stage
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).AMPLITUDE =zeros(size(AMPLITUDE(:)')) ;  % % Translation (centroid face)
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
DIRICHLET(icond).PRESCRIBED_DISP(iloadstate).INTERVAL =  [t0,tEND];  %
end