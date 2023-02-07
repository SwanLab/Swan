%% To-do
% (0) Mesh refining P1Function (RemeshingTests)
%       - Looks like mF is NOT exactly continuous...
%       - Perhaps remesh inside the refine method in P1Function?
% (1) Macro printing with FeFunctions
% (2) Micro printing with FeFunctions
% (3) Micro TopOpt printing

% EXTRAS
%  - Investigate: converting data to binary format to save read'n'write
%                 resources for paraview
%  - Tutorial for printing
%  - Study file ouptut size vs time (GiD/Paraview) to see which is better
%    for printing (test + graph)
%  - Check XY component of fgaussfunctions

%% Mid-term
% - Use FeFunctions in TopOpt_Problem
% - PDE belongs to Optimizer, not ShapeFunctional
% - Micro as three elasticity problems

%% Long-term
% - Move Input folder to a separat repository
% - Geometry only in Mesh
% - kill Mesh_Total (UnfittedMesh). still used somewhere but should be
%   removed
% - Recuperar gid unfitted mesh photo GiDimagecapturer
%      density --(project)--> unfittedmesh -> innermesh/photo

% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file