%% To-do
% DISCONTINUOUS MESH
% (0) Check projectors
% (1) Delete discontinuous meshes
%       IMPOSSIBLE
%           - P1DiscontinuousFunction
%           - Mesh
%           - Field
%           - MinimumDiscGradFieldWithVectorInL2
%               - used for fields
%       TOUGH
%           - DehomogenizingSingularitiesTest
%           - RemeshingTests
%           - Remesher -> make it disappear from everywhere? only in mesh?
%           ? SymmetricContMapCondition
%               - used for connecs, coords
%               x CoherentOrientationSelector
%           x LevelSetPeriodicAndOriented
%       EASY
%           x SingularitiesFinder

% (2) Delete fields
% (3) Delete Mesh_Total

%% Results
% FeFunction.project(target)

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

%% Backlog
% EXTRAS
%  - Investigate: converting data to binary format to save read'n'write
%                 resources for paraview
%  - Tutorial for printing
%  - Study file ouptut size vs time (GiD/Paraview) to see which is better
%    for printing (test + graph)
%  - Check XY component of fgaussfunctions