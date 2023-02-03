%% To-do
% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file



%%

% fem.print
%   - use functions
% topopt.print
%   - use functions
%   - optimizer: designVar, pde, CC (which have functionals)
%           - designvar pde and functionals should print
%           - create fefunctions in each of these classes for printing
%           *only*
%           - later on we will eventually do more stuff


% 
% - check XY componennt of fgaussfunctions
% - study file output size vs time (gid/paraview) to see which is better
%   for printing
%       - instructions at the wiki (running the test + graph)
% - update to master after cleanup

% - {{no fet}} uml de tot plegat (es una classe)
% - recuperar gid unfitted mesh photo GiDimagecapturer
%      density --(project)--> unfittedmesh -> innermesh/photo

% opengl, llicencies, gid 
% https://www.youtube.com/watch?v=zikDxtlvbUA

% ppt -> add gifs, not pictures

%% Comments

% Cantileverbeam_Hexahedra_Linear_Structured_SYM
% refining in P1DiscontinuousFunction vs P1Function (see RemeshingTests)

%% per mes endavant
% - moure input a repositori a part per alleugerir
% - geometry only in mesh
% - kill Mesh_Total (UnfittedMesh). still used somewhere but should be
%   removed