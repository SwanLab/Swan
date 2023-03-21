%% To-do
% DISCONTINUOUS MESH
% - Fix P1Function and P1Discontinuous function computeGradient fgauss
% nelem order
% - Should Fgaussdiscontfunctions exist on their own? RHSintegrator -> how
% do you assemble using this type of function? it just comes from the
% gradient of a p1fun / p1dfun

%% Questions
% - Filter_P1_LevelSet.getP0fromP1() ??
% - Careful: unfittedmesh needed to properly integrate!!!

%% Results
% 

%% Mid-term
% - Use FeFunctions in TopOpt_Problem
% - PDE belongs to Optimizer, not ShapeFunctional
% - Micro as three elasticity problems

%% Long-term
% - Move Input folder to a separate repository
% - Geometry only in Mesh
% - Recuperar gid unfitted mesh photo GiDimagecapturer
%      density --(project)--> unfittedmesh -> innermesh/photo

% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file

%% Backlog
% EXTRAS
%  - Investigate: converting data to binary format to save read'n'write
%                 resources for paraview
%  - Tutorial for printing
%  - Study file ouptut size vs time (GiD/Paraview) to see which is better
%    for printing (test + graph)
%  - Check XY component of fgaussfunctions