function [DATAOUT] = Q8_for_beam_18tests(AMPLITUDE_alltests,DATAOUT,ALL_surfaces_index,index_SURFACES_END_BEAMS)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% JAHO, 10-Apr-2024, Honest Green, C/Tuset, Barcelona
% See archetypal training script:
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE/Train_Q8_BEAMbeh.m
%------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end
 
%
% We begin by training with prescribed displacements
% We mimic the steps given in
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/Q8_for_beamBEHAVIOR_trainingTESTS.m

INPUTS_PARAMETERS = [] ;
INFO_BASIC_DEFORMATIONAL_MODES = [] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THESE ARE THE 6 TESTS WITH PRESCRIBED RIGID BODY DISPLACEMENT ON ONE OF
% THE ENDS
% --------------------------------------------------------------------------------------------------
AMPLITUDE = AMPLITUDE_alltests*eye(7) ;
AMPLITUDE = AMPLITUDE(:,1:6);
for itest = 1:size(AMPLITUDE,2)
    INPUTS_PARAMETERS(itest).DIRICHLET =  DirichBound_BEAM_generic(DATAOUT.t0,DATAOUT.tEND,AMPLITUDE(:,itest),index_SURFACES_END_BEAMS)  ;
    INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = 'RIGID_BODY_WARPING_ON_SURFACES'; %
    INPUTS_PARAMETERS(itest).NEUMANN = [] ;
end

% TESTS WITH PRESCRIBED FORCES ON both ends
% ------------------------------------------

AMPLITUDE_FORCE= 1e3*eye(6);
itest  = length(INPUTS_PARAMETERS) ;



ifixed = [1,2 ];
iforce = [2,1];
for isurfLOC = 1:length(ifixed)
    for itestLOC = 1:size(AMPLITUDE_FORCE,2)
        itest = itest + 1 ;
        % DIRICHLET CONDITION
        % ----------------------------------------------------------------------
        INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = 'RIGID_BODY_WARPING_ON_SURFACES'; %
        INPUTS_PARAMETERS(itest).DIRICHLET.PRESCRIBED_DISP.AMPLITUDE = zeros(1,7) ;
        INPUTS_PARAMETERS(itest).DIRICHLET.PRESCRIBED_DISP.TIMEFUN = @(t) t/DATAOUT.tEND ;  %
        INPUTS_PARAMETERS(itest).DIRICHLET.PRESCRIBED_DISP.INTERVAL =[DATAOUT.t0,DATAOUT.tEND];  %
        INPUTS_PARAMETERS(itest).DIRICHLET.NUMBER_SURFACE = index_SURFACES_END_BEAMS(ifixed(isurfLOC)); %
        % NEUMAN CONDITIONS
        INPUTS_PARAMETERS(itest).NEUMANN.NUMBER_SURFACE =index_SURFACES_END_BEAMS(iforce(isurfLOC)) ;  % Surface on which the load is applied
        INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE = [] ;
        INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.AMPLITUDE = AMPLITUDE_FORCE(:,itestLOC); % Force per unit surface
        INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.TIMEFUN =    @(t) t/DATAOUT.tEND ;  %
        INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.INTERVAL = [DATAOUT.t0,DATAOUT.tEND];  %
        INPUTS_PARAMETERS(itest).NEUMANN.GENERALIZED_POINT_FORCE.ISLOCAL = 0 ;
        
    end
end

DATAOUT.INPUTS_PARAMETERS = INPUTS_PARAMETERS;

INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_DEFORMATIONAL_MODES = 18 ;  % 24-6
DATAOUT.INFO_BASIC_DEFORMATIONAL_MODES = INFO_BASIC_DEFORMATIONAL_MODES;
DATAOUT.METHOD_TRAINING_ELASTIC_RANGE = 'INVARIANT_MODES';


% 
% 
% %
% %
% %
% % ini  = 0 ;
% 
% %
% %  for itestLOC = 1:size(AMPLITUDE,2)
% %     itest = ini + itestLOC ;
% %     INPUTS_PARAMETERS(itest).DIRICHLET.AMPLITUDE = AMPLITUDE(:,itestLOC) ;
% %     INPUTS_PARAMETERS(itest).DIRICHLET.TIMEFUN = @(t) t/DATAOUT.tEND ;  %
% %     INPUTS_PARAMETERS(itest).DIRICHLET.INTERVAL =[DATAOUT.t0,DATAOUT.tEND];  %
% %     INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = 'FE_SHAPE_functions_only_selected_surfaces'; %
% %     INPUTS_PARAMETERS(itest).DIRICHLET.BOUNDARY_SURFACES_TO_INCLUDE = index_SURFACES_END_BEAMS; %
% % end
% 
%  ntests = 24;
%  AMPLITUDE = AMPLITUDE_alltests*eye(ntests) ;
% ini  = length(INPUTS_PARAMETERS) ;


%INFO_BASIC_DEFORMATIONAL_MODES.INDICES_for_invariant_TRAINING_TESTS= 1:ini;  %length(INPUTS_PARAMETERS) ;

% 
% for itestLOC = 1:size(AMPLITUDE,2)
%     itest = ini + itestLOC ;
%     INPUTS_PARAMETERS(itest).DIRICHLET.AMPLITUDE = AMPLITUDE(:,itestLOC) ;
%     INPUTS_PARAMETERS(itest).DIRICHLET.TIMEFUN = @(t) t/DATAOUT.tEND ;  %
%     INPUTS_PARAMETERS(itest).DIRICHLET.INTERVAL =[DATAOUT.t0,DATAOUT.tEND];  %
%     INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = 'FE_SHAPE_functions_only_selected_surfaces'; %
%     INPUTS_PARAMETERS(itest).DIRICHLET.BOUNDARY_SURFACES_TO_INCLUDE = ALL_surfaces_index; %
%     INPUTS_PARAMETERS(itest).NEUMANN = [] ;
%     
% end
% 
% 
% INFO_BASIC_DEFORMATIONAL_MODES.INDICES_all_prescribed_TRAINING_TESTS= (ini+1):length(INPUTS_PARAMETERS) ;   %length(INPUTS_PARAMETERS) ;
% 


