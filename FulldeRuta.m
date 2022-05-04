%% To-do
% a) DIMENSIONS
%       OK! - Delete nunkn in RHS
%       OK! - Avoid using nelem from dim
%       OK! - Avoid using ngaus from dim 
%       OK* - nstre only in elasticity  (maybe in Bmatrix) and rename it
%             as nVoigt
%       IDK - ndofPerNode (and ndofPerElem) in interpolation times npnod
%             (nnode) of mesh gives ndof of field (in dim)
%                   - ndofPerNode in interpolation?
%       OK* - nnode ---> nnodePerElem (in Mesh)
%       OK! - npnod ---> nnodes (in Mesh)
%       OK! - create dimensions from s.type (Scalar/Vector)

% b) INTEGRATORS
%       OK* - "integrate" all RHS
%               - note: pending more cleanup.

% c) TESTS
%       NO* - Quadratic shape functions for thermal, elastic, elastic_micro
%             with corresponding tests
%               - When the interpolation is set as quadratic, the
%                 connectivity matrix must be once again calculated.
%                 ConnecCoordFromInterpAndMesh is supposed to do it, but it
%                 does not work as intended.
%               - Some elements overlap, the original node numeration is
%                 lost in the process and lose physical meaning
%               - As a result (?), the stiffness matrix is not invertible.
%                 If the dofConnec is not to be calculated, the
%                 alternatives may be more expensive.

%         pd.bc.dirichlet = [10 1 0; 10 2 0; 7 1 0; 7 2 0; 19 1 0; 19 2 0];
%         pd.bc.pointload = [24 2 -1];

% d) STOKES
%       WIP - Restore Stokes_Problem
%               - FemTestsSuite, FemTests
%               - StokesComputer, StokesFEM

% z) LONG-TERM
%       WIP - Re-use FeFunction for displacements... with its own
%             dimensions

% delete problemData
% refactoring DImensions (mesh in init and compute)
% delete nnodes, nnodeElem from Dimensions (or private)