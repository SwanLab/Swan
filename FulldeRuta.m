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

