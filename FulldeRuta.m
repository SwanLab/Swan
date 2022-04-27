%% To-do
% a) DIFFREACTPROBLEM
%       OK! - In Filter use LHSintegrator rather than DiffReact.
%       OK! - Delete getM and getK and setEpsilon and computeDvolume
%       OK* - setLHStype in DiffReact in Filter
%       OK! - DiffReactTestsSuite created

% b) ELASTICPROBLEM, DIMENSIONS and INTEGRATORS
%       OK! - delete interp in ElasticProblem
%       OK* - dimEscalar, vector....by fields...nElem,nDim private in dim
%               - DimensionScalar, DimensionVector
%       OK* - Create RHSintegrator for Elastic, ElasticMicro and thermal

        % Comments:
        %       - I started inverting the approach to Integrators and
        %         RHSintegrators, as discussed a while ago. Previously,
        %         RHSintegrators computed only the elemental RHS via
        %         fGauss, and *not* the RHS.
        %       - Is it useful, though?
        %       - MinimumGradFieldWithVectorInL2, called from
        %         DilationFieldComputer. fNodal?
        %       - Tests are missing to cover stuff eg. HarmonicProjector.
        %         HarmonicProjectionExample exists, but...

% d) STOKES
%       WIP - Restore Stokes_Problem
%               - FemTestsSuite, FemTests
%               - StokesComputer, StokesFEM


% fix newdimensions in ElasticProblem
% quadrature as input in strain
% in Filter create MassMatrix when construction


% +++QUadratic shape functions for thermal, elastic, elastic_micro with
% corresponding tests
% ++++"integrate" all RHS
% ++++nstre only in elasticity  (maybe in Bmatrix) and rename it as nVoigt
% +++trying to avoid using nelem from dim
% +++trying to avoid ngaus from dim 
% ++++ndofPerNode (and ndofPerElem) in interpolation times npnod (nnode) of mesh gives ndof of field (in
% dim)
% ++++ nnode ---> nnodePerElem (in mesh)
% create dimensions from type
% delete nunkn in RHS

% Re-use FeFunction for displacements.... with its own dimensions
