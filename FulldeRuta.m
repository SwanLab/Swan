%% To-do
% a) DIFFREACTPROBLEM
%       OK! - merge diffreacProblemRobin and Neuman...
%       OK!      - dispacth trhough LHSintegrator
%       OK!      - merge problemdimensions and dimensions, no need
%       OK! - create LHS integrator for boundaryMassmatrix
%       YET - DiffReactProblem tests

        % Comments:
        %       - DiffReactTests: RHS? Dirichlet? Neumann?
        %       - Separate problemLHS from LHSintegrators in a new class?

% b) MINOR CLEANUP
%       YET - improve performance graph using average + deviation
%       OK! - Precomputedvariabletest overwrite results

        % Comments:
        %       - AbstractSettings was the one, line 89

% c) ELASTICPROBLEM
%       OK! - masterSlaveNodes not in Mesh...in BoundaryConditions
%                - Pending createRectangularMesh at ShapesInMicrostructures
%       OK! - s.material as an input (FemDataContainer)
%       YET - refactoring ElasticProblemMicro

        % Comments:
        %       - Perhaps material should not take ngaus
        %       - This way, quad, interp, and geometry *may* only be used
        %         as part of the RHS + strain/stress
        %       - If we end up doing the problemLHS thing, it could also
        %         lead to a problemRHS with all of this

% d) STOKES
%       WIP - Restore Stokes_Problem
%               - FemTestsSuite, FemTests
%               - StokesComputer, StokesFEM

        % Comments:
        %       - Seems to be made just for that one test. See
        %         test2d_stokes_triangle (Vol_force, pressure, velocity)
        %       - Is the mesh even used? Looks like it's all interpolations
        %       - Separate DimensionVariables? Separate BoundaryConditions?
        %       - The RHS thing could be useful here. The ForcesComputer in
        %         ElasticProblem should be completely different, maybe a
        %         change in the architecture could make it easier.

% e) PREPROCESS
%           - Refactoring PreProcess -- FemInputReader_GiD.m



% OK! In Filter use LHSintegrator rather than DiffReact
% OK! DiffReact delete getM and getK and setEpsilon and computeDvolume
% OK* setLHStype in DiffReact in Filter
% YET Create RHSintegrator for Elastic, ElasticMicro and thermal
% OK! delete interp in ElasticProblem
% YET dimEscalar, vector....by fields...nElem,nDim private in dim
% YET Stokes..

% Check test_VigergauzMicroStructureFromStrain, test_bridge
