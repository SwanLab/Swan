function [NEUMANN] = Neumann_OPPrb_generic1(t0,tEND,AMPLITUDE,SURFACES)
% 
% 
%   INPUTS_PARAMETERS(itest).NEUMANN.NUMBER_SURFACE =index_SURFACES_END_BEAMS(iforce(isurfLOC)) ;  % Surface on which the load is applied
%         INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE = [] ;
%         INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.AMPLITUDE = AMPLITUDE_FORCE(:,itestLOC); % Force per unit surface
%         INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.TIMEFUN =    @(t) t/DATAOUT.tEND ;  %
%         INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.INTERVAL = [DATAOUT.t0,DATAOUT.tEND];  %
%         INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.ISLOCAL = 0 ;


 
%%%%%%

icond = 1; % Number of condition
%---------------------------------------------------------------
NEUMANN(icond).NUMBER_SURFACE = SURFACES(1) ;   % Number of SURFACE on which DISPLACEMENT  is prescribed
iloadstate = 1;  % Loading stage
NEUMANN(icond).GENERALIZED_POINT_FORCE(iloadstate).AMPLITUDE =AMPLITUDE(:)' ;  % % Translation (centroid face)
NEUMANN(icond).GENERALIZED_POINT_FORCE(iloadstate).TIMEFUN =  @(t) t/tEND ;  % Function definining the temporal evolution
NEUMANN(icond).GENERALIZED_POINT_FORCE(iloadstate).INTERVAL =  [t0,tEND];  %
