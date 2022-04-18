%% To-do
% a) DIFF-REACT
%       OK! - merge diffreacProblemRobin and Neuman...
%       OK!      - dispacth trhough LHSintegrator
%       OK!      - merge problemdimensions and dimensions, no need
%       OK! - create LHS integrator for boundaryMassmatrix
%       YET - DiffReactProblem tests
%       YET - Maybe some additional cleanup can be made

% b) MINOR CLEANUP
%       YET - improve performance graph using average + deviation
%       OK! - Precomputedvariabletest overwrite results

% c) ELASTIC
%       OK! - masterSlaveNodes not in Mesh...in BoundaryConditions
%                - Pending createRectangularMesh at ShapesInMicrostructures
%       OK! - s.material as an input (FemDataContainer)
%       YET      - Perhaps material should not take ngaus
%       YET - refactoring ElasticProblemMicro

% d) STOKES
%       YET - Restore Stokes_Problem
%               - FemTestsSuite, FemTests
%               - StokesComputer, StokesFEM

% e) PREPROCESS
%           - Refactoring PreProcessing--FemInputReader


% masterSlaveNodes not in Mesh...in BoundaryConditions
% diffreact tests
% create LHS integrator for boundaryMassmatrix
% Precomputedvariabletest overwrite results
% merge diffreacProblemRobin and Neuman...dispacth trhough LHSintegrator
% createMaterial
% Refactoring ElasticProblemMicro
% Restore Stokes
% Refactroing PreProcessing--FemInputReader