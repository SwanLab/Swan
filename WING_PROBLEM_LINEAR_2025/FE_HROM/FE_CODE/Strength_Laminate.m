function   DATAOUT =   Strength_Laminate(DATAOUT,MATERIAL,DATA)
%dbstop('3')
if nargin == 0
    load('tmp.mat')
end
% --------------------------------------------------------------
% Computing  ULTIMATE BENDING AND TWIST MOMENTS   FOR LAMINATES
% --------------------------------------------------------------

% COmpute average element stress for the elements of each ply
% ----------------------------------------------------
stressPLY = StressElementAvg_ply(DATAOUT) ;

% Number of plies
%-----------------
nplies = length(unique(DATAOUT.MaterialType)) ;
%
MaximumValueFailure = zeros(nplies,1) ;
StressRotated = {} ;
for iply = 1:nplies
    % --------------------------------------------------------------------
    % For each ply, we calculate the stresses in the ply referency system
    % --------------------------------------------------------------------
    % Stress rotation matrix
    angle = MATERIAL.PLY(iply).ANGLE ;
    T = RotationMatrix(angle) ;
    % We rotate all stresses
    stresses = inv(T)*stressPLY{iply};
    
    % We recover the matrix defining the failure envelope of the composite
    % material
    load(MATERIAL.PLY(iply).NAMEWS_FAIL) ;
    % Failure criterion considering only membrane stresses
    comp = [1 2 6] ;
    stresses = stresses(comp,:) ;
    F_stress = Fenv*stresses ;
    FailureCrit = sum(stresses.*F_stress,1) ;
    nfail = length(FailureCrit) ;
    nfail =  min(nfail,DATA.NUMBER_ELEMENTS_FAIL);% failed
    % Sort -->
    FailureCrit = sort(FailureCrit,'descend') ;
    % Maximum value of stres'*F*stress
    MaximumValueFailure(iply) = sqrt(FailureCrit(nfail));
    
    DATA= DefaultField(DATA,'PLOT_FAILURE_POINTS',0) ; 
    if DATA.PLOT_FAILURE_POINTS == 1
     PlotStressPoints(iply,Fenv,stresses) ;

        %
        
    end
    
    
end


% The   load factor to be employed is the maximum of MaximumValueFailure
% ------------------------------------------------------------------------
[load_factor iplyFAIL]= max(MaximumValueFailure) ;

% Finally, hence, the valu of the MACRO-stress vector when the RVE starts to fail
% is
stressAVG_FAIL = DATAOUT.stressAVG/load_factor ;
% Likewise, for strains

DATAOUT.stressAVG_FAIL = stressAVG_FAIL ;
DATAOUT.iplyFAIL = iplyFAIL ;
DATAOUT.LOAD_FACTOR_FAIL = load_factor ;
%DATAOUT.strainAVG_FAIL = strainAVG_FAIL ;
