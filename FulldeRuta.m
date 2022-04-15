%% To-do
% a) CLEANUP
%       OK! - Rename New stuff to just the name
%                - ElasticProblem, DiffReactProblem
%                - ElasticProblemMicro, DiffReactProblemMicro
%       OK! - Delete old DOFs and Elements
%                - Stokes, Hyperelastic?
%       OK! - Move dofsInElem to LHSintegrator
%           - Stokes_Problem still uses BoundaryConditionsApplier
%                - Left for legacy purposes

% b) CLOSING
%       OK! - Adapt ThermalProblem via NewDiffReactProblem
%       OK! - Add vars2print for micro

% c) RESULT VISUALIZATIONS
%       OK! - Set up GiD for ThermalProblem and visualize results
%       OK! - Set up GiD for Micro problem and visualize results

% d) PERFORMANCE
%           - See report.mlx


%%% 
% OK! replace applyNode --- for nnode = size(globalconnec,2) in Assembler
% OK! Clean dofsInElem THerma, Elastic
% OK! dofInElem by composition --->  Assembler and Forces
% OK! MaterProperties data as input 
% OK! Performance (time) vs nElem: for different cases

% Clean ElasticProblemMicro ---> vars, vars2Print, Chomog, variables....

% OK! DiffReac in Newmann and Robin (Two diff react)

% Clean Thermal, Elastic, DiffReac --- (only) dim, boundary, LHS,RHS, u..
%   - DiffReactProblem: it's pretty clean, since quad, interp and geom are
%                       not called anywhere else.
%   - ElasticProblem: it's tough, since it really needs to be calculated in
%                     that class as it's needed elsewhere. Moved to init.


% Test coming from GiD (Swan.gid)

% masterSlaveNodes not in Mesh...in BoundaryConditions
% diffreact tests
% create LHS integrator for boundaryMassmatrix
% Pcomputed variable rewritte results
% merge diffreacProblemRobin and Neuman...dispacth trhough LHSintegrator
% createMaterial
% Refactoring ElasticProblemMicro
% Refactroing PreProcessing--FemInputReader