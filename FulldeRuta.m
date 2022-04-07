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
%dofInElem by composition --->  Assembler and Forces
% DiffReac in Newmann and Robin (Two diff react)
% Clean Therma, Elastic, DiffReac --- (only) dim, boundary, LHS,RHS, u..
% Clean dofsInElem THerma, Elastic
% MaterProperties data as input 
% Clean ElasticProblemMicro ---> vars, vars2Print, Chomog, variables....
% Test coming from GiD (Swan.gid)
% Performance (time) vs nElem: for different cases
% replace applyNode --- for nnode = size(globalconnec,2) in Assembler


