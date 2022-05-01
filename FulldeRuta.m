%% To-do
% a) DIMENSIONS
%       OK! - Delete nunkn in RHS
%       OK! - Avoid using nelem from dim
%       OK! - Avoid using ngaus from dim 
%       OK* - nstre only in elasticity  (maybe in Bmatrix) and rename it
%             as nVoigt
%       YET - ndofPerNode (and ndofPerElem) in interpolation times npnod
%             (nnode) of mesh gives ndof of field (in dim)
%       YET - nnode ---> nnodePerElem (in Mesh)
%       YET - create dimensions from s.type (Scalar/Vector)

% b) INTEGRATORS
%       YET - "integrate" all RHS

% c) TESTS
%       YET - Quadratic shape functions for thermal, elastic, elastic_micro with
%             corresponding tests

% d) STOKES
%       WIP - Restore Stokes_Problem
%               - FemTestsSuite, FemTests
%               - StokesComputer, StokesFEM

% z) LONG-TERM
%       WIP - Re-use FeFunction for displacements... with its own
%             dimensions


%     - Quadratic shape functions for thermal, elastic, elastic_micro with
%       corresponding tests
%     - "integrate" all RHS
% OK* - nstre only in elasticity  (maybe in Bmatrix) and rename it as nVoigt
% OK! - trying to avoid using nelem from dim
% OK! - trying to avoid ngaus from dim 
%     - ndofPerNode (and ndofPerElem) in interpolation times npnod (nnode)
%       of mesh gives ndof of field (in dim)
%     - nnode ---> nnodePerElem (in Mesh)
%     - create dimensions from s.type (Scalar/Vector)
% OK! - delete nunkn in RHS

% Re-use FeFunction for displacements.... with its own dimensions
