function [INPUTS_PARAMETERS] = Q8_for_beamBEHAVIOR_ALL(AMPLITUDE_alltests,DATAOUT)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
% JAHO, 26-March-2024, plane from Barcelona to Tenerife
%------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end

INPUTS_PARAMETERS = [] ;
INFO_BASIC_DEFORMATIONAL_MODES = [] ; 

nDOFS = 24; 
AMPLITUDE = eye(nDOFS) ;
for itest = 1:size(AMPLITUDE,2)
    INPUTS_PARAMETERS(itest).DIRICHLET.AMPLITUDE = AMPLITUDE(:,itest)        ;
    INPUTS_PARAMETERS(itest).DIRICHLET.TIMEFUN = @(t) t/DATAOUT.tEND         ;  %
    INPUTS_PARAMETERS(itest).DIRICHLET.INTERVAL =[DATAOUT.t0,DATAOUT.tEND]   ;  %
    INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = ''           ; %
end



INFO_BASIC_DEFORMATIONAL_MODES.NUMBER_OF_BASIC_DEFORMATIONAL_MODES = length(INPUTS_PARAMETERS) ; 

disp(['Q8 EIF ELEMENT FOR  representing beam behavior'])
disp(['TOTAL NUMBER OF TESTS = ',num2str(length(INPUTS_PARAMETERS) ),' = ', 'number of deformational modes'])


% else
%     %---------------------------------------------------------------
%     AMPLITUDE = eye(24)*AMPLITUDE_alltests ;
%     ini  = length( INPUTS_PARAMETERS) ;
%     for itestLOC = 1:length(AMPLITUDE)
%         itest = ini + itestLOC ;
%         INPUTS_PARAMETERS(itest).DIRICHLET.AMPLITUDE = AMPLITUDE(:,itestLOC) ;
%         INPUTS_PARAMETERS(itest).DIRICHLET.TIMEFUN = @(t) t/DATAOUT.tEND ;  %
%         INPUTS_PARAMETERS(itest).DIRICHLET.INTERVAL =[DATAOUT.t0,DATAOUT.tEND];  %
%         INPUTS_PARAMETERS(itest).TypeFunctionDisplacementTRAINING = ''; %
%     end
%
% end

