classdef PhaseFieldTests < handle & matlab.unittest.TestCase

% En el futur: agrupar testsTO malla en comú i sense definir functionals, opt..;
% crear difs funcions a topopttests depenent tipus functionals...

    properties (TestParameter)
        testsTO = {
            'test_cantilever', 'test_cantilever2', 'test_cantilever3', ...
            'testDualNestedInPrimal_WithProjectedGradient', ...
            'testDualNestedInPrimal_WithSlerp', ...
            'test_interiorPerimeter',...
            'test_anisotropy','test_anisotropy_interior','test_nullspace',...
            'test_interiorPerimeterPDErho','test_filterLump','test_cantilever_IPM',...
            'test_dirichletProjection','test_gripping','test_micro', 'test_micro2'
            }
    end

    methods (Test, TestTags = {'PF'})

        function test1DAmbrosioTortorelli(testCase, testsTO)
            run(testsTO);
            s.

            problem = TestingPhaseField(s);
            ..

            phiNew =
            load([testsTO,'.mat'],'x');
            err = norm(x - xNew)/norm(x);
            tol = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol);
        end

    end
end