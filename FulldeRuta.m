%% To-do
% a) THE GOOD

%       OK! - delete problemData in ElasticProblem
%       OK! - refactoring Dimensions (mesh in init and compute)
%       MEH - delete nnodes, nnodeElem from Dimensions (or private)

%           - NewStokesProblem now looks a bit better. Changes made to the
%             input (via StokesDataContainer) and the solver (via
%             NonLinear_Solver)
%           - ConnecCoordFromInterpAndMesh appears to work mostly* fine. It
%             is completely fine for sure in the Stokes case (hand-checked)

% b) THE BAD

%           - *For the simpleInterpTest, the globalConnec given is wrong,
%             but it may be due to it being a two-element mesh.
%           - For quadratic interpolations, the stiffness matrix is
%             singular, and I have no clue as to what is going on
%             input (via StokesDataContainer) and the solver (via
%             NonLinear_Solver)
%                   - Created InterpolationTranslator to solve the mismatch
%                     between pre-interpolation and post-interp data.
%                   - The assembly seems to have been performed properly
%                   - The boundary conditions seem to have been applied
%                     correctly. Adding an additional Dirichlet condition
%                     does not seem to help.

% c) THE UGLY

%           - ndof and nnode is still needed in some cases, mainly due to
%             issues with ndof and Boundary Mass Matrices. Ambitious fix
%             proposed below
%           - Physical meaning of the matrices at Element_Stokes? Move to
%             LHSintegrator?
%           - What to do with abandoned methods at LHSintegrator + old
%             assembly?

% Possible concept?
% Instead of including DimensionVariables as part of the problem, make a
% new class Field which is part of the problem. THIS class could have the
% dimensions of the field, its own interpolation + connectivities, and
% eventually, its own FeFunction

% d) STOKES
%       WIP - Restore Stokes_Problem
%               - FemTestsSuite, FemTests
%               - StokesComputer, StokesFEM

% z) LONG-TERM
%       WIP - Re-use FeFunction for displacements... with its own
%             dimensions

%- Tests quadratic
%- Stokes