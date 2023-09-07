%% To-do

% WEBSITE
    % Website: think of something to easily add new simulations
    % Readme: take out everything, replace with website, including some
    %         examples, include link to website
    
% TUTORIALS
    % - 0. Create mesh (Ton)
    % - 1. Fem Thermal (Ton) (include plots & prints)
    % - 2. Fem Elasticity (Ton) (include plots & prints)
    % - 3. Create level set functions (square with circle) (Ton)
    %      Extrude+export stl tutorial
    % - 4. Filter (Jose)
    % - 5. TopOpt (Jose)
    %   5.1. Macro
    %   5.2. Micro
    % - 6. Shape optimization (Alex)
    % - 7. Dehomogenization (Alex)

% EXTRAS
    % - Projectors (Jose)
    % - Subdomains (Lagrange Multipliers)
    % - Micro elasticity
    % - Optimizers

% Nesterov? + Arnau


%% Comments

%% Questions

%% Results


%% Long-terms
% Filters -> using only LHS/RHsinteg
% RHS with test and trial

% CharacteristicFunction should return an UnfittedMesh

% Nesterov


% Preprocess / GUI -> TBD


% - Use FeFunctions in TopOpt_Problem

%% Backlog
% - Move Input folder to a separate repository
% - Geometry only in Mesh

% EXTRAS
%  - Investigate: converting data to binary format to save read'n'write
%                 resources for paraview
%  - Study file ouptut size vs time (GiD/Paraview) to see which is better
%    for printing (test + graph)
% - Recuperar gid unfitted mesh photo GiDimagecapturer
%      density --(project)--> unfittedmesh -> innermesh/photo

% - GiDImageCapturer -> some variables should be "user variables", defined
%   once in a centralized file