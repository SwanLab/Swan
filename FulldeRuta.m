%% Updates
% - Stokes: For some reason, Transient is consistently way faster than
%           Steady.
% - Field: BoundaryConditions no longer there, but they can be used to
%          translate.
% - LHSintegrators: they now take Fields as inputs.


%% To-do
% - Move BC from fields
%       - many possibilities to do it, but none that don't smell.
%           - BoundaryConditionFactory that generates BoundaryCondition,
%             just like filters?
%           - Passing a cell array as the input BC for the Boundary
%             Conditions class? (you lose info like ndimfield)
%           - Passing a struct where each field indicates type/domain/...?
%       - It makes no sense to have free_dofs for each field in a cell
%         array. BoundaryConditions were moved away from Field because of
%         that.
%       - it is difficult to sense which will be the winner, perhaps it's
%         time to push with new features and see what comes up.
% - move geometry to Mesh
%       - it cannot work properly. For problems such as stokes, there are
%       two different interpolations. creating the geometry inside the mesh
%       only allows for one interpolation at a time -- and you cannot copy
%       the mesh to make a new one becaues of the way matlab handles it
% - give importance to interptranslator in field
% - LHSintegrator types: field/test function

%% Long-term
% - FeFunction and Field should converge into one. Mesh.coord could be a
%   FeFunction: given a quadrature, it returns xGauss 

% - Try out ways to make constructors more elegant, perhaps by grouping up
%   properties and looping over them, some tagged as RequiredProps and the
%   rest as OptionalProps.

%% Previous To-do

%       - The new way of computing the stiffness matrix is
%         significantly faster (previously: 0.4s, now <0.1s). Also, the
%         material is basically irrelevant (identity matrix...). Now
%         uses the Assembler class (+ sym grad?)

% Field: - Is it overreaching?
%        - does a field have BCs? Even Neumann? (maybe, but if there are
%          two fields who applies them properly?)
%        - the need to rethink and refactor BCs comes up

% LHSintegrators using only Fields? Does it make sense physically in all
% cases (eg. ScalarProduct, ShapeFunctional, Poperator)? Also, Assembler...

% Element_Stokes.compute_D -> interp.ndime? Why pressure quadratic
% quadrature? Also, Matlab feature/issue with copying objects

% Any elegant way to avoid using if isfield to simplify constructors?

%% Previous
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