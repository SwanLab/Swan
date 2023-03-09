%% To-do
% DISCONTINUOUS MESH
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

% (2) Delete fields
%       - Note: ignored P2 elastic for now
%       - ELASTIC STIFFNESS
%           Pending test_anisotropy cleanup
%       - STIFFNESS
%           Pending: test1DLHS -> Geometry_Line to include computeInvJac
%           a la Geometry_Volumetric
%           Pending: diffreact dependencies
%       - MASS
%           Pending: test1DLHS, LHSintegrator_Stokes
%           Noteworthy: LHSintegrator_MassBoundary

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