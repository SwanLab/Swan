%% To-do
% (0) Mesh refining P1Function (RemeshingTests)
%       - Looks like mF is NOT exactly continuous...
%       - Perhaps remesh inside the refine method in P1Function?
% (DONE) Macro printing with FeFunctions
% (DONE) Micro printing with FeFunctions
%       - How do we print the three primal cases without three elastic
%         problems
% (DONE) Micro TopOpt printing
%       - Memory leak when printing LevelSet instead of Density...
%       - alpha becomes 0x1x0, perimeterp0 is a 6400x1x6400
%       - different output formats in computeVolumeFraction !!!

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
% - Move Input folder to a separate repository
% - Geometry only in Mesh
% - kill Mesh_Total (UnfittedMesh). still used somewhere but should be
%   removed
% - Recuperar gid unfitted mesh photo GiDimagecapturer
%      density --(project)--> unfittedmesh -> innermesh/photo

% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file