%% To-do
% - Check Paraview micro level-set
% - Check STL 3D
% - mass p1p0

%% Questions

%% Results
% 


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