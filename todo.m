%% To-do
    % Fisic: open hole // grans vs petites. trobar llibre document
    % referencia
    % Convergenece: monolithic -> reduced
    % No voigt: matrius geometrica i no lineal (mathcha + implementar)
    % Stretches
    % Open hole: buscar exemple deformada sencera, no nomes un quart


    % Speed-up: pass ndofs to BCApplier, calculate free dofs just once
    
% WEBSITE
    % STEP 1
    % - Readme: tot fora, portar cap a la web + documentation
    %           afegir nous exemples visuals per cridar atencio
    % - Web: nous exemples, afegir manuals/tutorials tipus Firedrake
    % - Fer tutorial TopOpt (Cantilever compliance -> gif + monitoring)
    %   -> començar a pensar amb settings

    % STEP 2
    % FEM nou seguint model idealWorld.m

    % STEP 3
    % TopOpt nou seguint model idealWorld.m

% TUTORIALS
    % - 0. Create mesh (Ton)
    % - 1. Fem Thermal (Ton) (include plots & prints)
    % - 2. Fem Elasticity (Ton) (include plots & prints)
    % - 3. Create level set functions (square with circle) (Ton)
    %      Extrude+export stl tutorial
    % - 4. Filter (Jose)
    % - 5. TopOpt 
    %   5.1. Macro (Ton)
    %   5.2. Micro (Jose)
    % - 6. Shape optimization (Alex)
    % - 7. Dehomogenization (Alex)
    % - 8. Projectors (Ton)

% EXTRAS
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