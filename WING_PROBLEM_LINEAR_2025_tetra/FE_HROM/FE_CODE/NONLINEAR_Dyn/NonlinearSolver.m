function NonlinearSolver(OPERfe,DATA,DynDATA,PROPMAT)
% INPUTS
% - OPERfe: FE operators, and sundry variables
% - DATA: Various inputs data.
% - DynDATA: Input data for dynamic problems
% PROPMAT : Initalization internal variables, internal parameters 
% JAHO, 17-Sept-2018. Copy of
%/home/joaquin/Desktop/CURRENT_TASKS/POTENTIAL_RESEARCH_TOPICS/COMBINING_MULTISCALE_REDUCTIONMODELS/REPORT_MULTIS_REDUC_MODEL/PAPER_02/MATLAB_FILES/FE2DAnalysisDyn.m

if  nargin == 0
    load('tmp.mat')
end

