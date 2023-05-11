%% To-do
% - Check Paraview micro level-set
% - Check STL 3D
% - mass p1p0

% 1) ElasticProblemMicro_Fast
%       - one loop at solve(). displacements+strain+stress+fluct calculated
%         without adding additional loops
%       - code cleanup (see computeDisplacements())
%       - delete vars/variables/whatever -> change functionals

% 2) FunctionPrinter/ParaviewPostprocessor -> move to Factory
% 3) Mesh.print

%% Questions

%% Results
% 

% .setC 
% OrientationAsDesignVariable
% ShFunc_Chomog
% ShFunc_Compliance
% ShFunc_ComplianceComparison_constraint
% ShFunc_Compliance_constraint
% ShFunc_NonSelfAdjoint_Compliance
% ShFunc_StressNorm
% ShFunc_StressNorm2
% ShFunc_StressNorm3
% testComputingFemWithVademecumData


%% Long-term
% Mesh cleanup public
% Trial/test -> Mass matrix P0 P1
% Filters -> using only LHS/RHsinteg

% CharacteristicFunction should return an UnfittedMesh

% Nesterov


% Preprocess / GUI -> TBD


% - Use FeFunctions in TopOpt_Problem
% - PDE belongs to Optimizer, not ShapeFunctional
% - Micro as three elasticity problems

%% Backlog
% - Move Input folder to a separate repository
% - Geometry only in Mesh

% EXTRAS
%  - Investigate: converting data to binary format to save read'n'write
%                 resources for paraview
%  - Tutorial for printing
%  - Study file ouptut size vs time (GiD/Paraview) to see which is better
%    for printing (test + graph)
%  - Check XY component of fgaussfunctions
% - Recuperar gid unfitted mesh photo GiDimagecapturer
%      density --(project)--> unfittedmesh -> innermesh/photo

% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file