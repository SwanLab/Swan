%% To-do
% OK! - Delete rhs_shapefunction -> use funs
% OK! - Clean RHSunfitted
% OK! - switch class(mesh.type) @ RHSfactory
% OK! - filter_pde_density using rhs_unfitted (now rhs_shapefunctionfun*
% OK! - filter_pde_levelset name change for unfitted

%% Questions
% - Filter_P1_LevelSet.getP0fromP1() ??
% - Careful: unfittedmesh needed to properly integrate!!!

%% Results
% 


%% Long-term
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