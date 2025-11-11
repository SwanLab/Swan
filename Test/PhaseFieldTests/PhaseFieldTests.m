classdef PhaseFieldTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        singleElement = {'AT1','AT2','Force','Split','HomogGrad'}
        singleElementForce = {'Force'}
        complexCases = {'SEN'}
        homogenizationCases = {'Square','Hexagon','Ellipse'}
    end

    methods (Test, TestTags = {'PF'})
        function testPhaseFieldSingleElem(testCase,singleElement)
            filename = ['testPhaseFieldSingleElem',singleElement];
            load(filename,'input');
            tester = TestingPhaseField(input);
            outputData = tester.compute();
            xNew.F = outputData.force;
            xNew.U = outputData.displacement.value;
            load(filename,'xRef');
            errF = norm(xNew.F-xRef.F)/norm(xRef.F);
            errU = norm(xNew.U-xRef.U)/norm(xRef.U);
            tol  = 1e-6;
            testCase.verifyLessThanOrEqual(max(errF,errU), tol)
        end

        function testPhaseFieldComplexCases(testCase,complexCases)
            filename = ['testPhaseField',complexCases];
            load(filename,'input');
            tester = TestingPhaseField(input);
            outputData = tester.compute();
            xNew = outputData.force;
            load(filename,'xRef');
            err = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testPhaseFieldHomogenization(testCase,homogenizationCases)
            filename = ['testPhaseFieldHomogenization',homogenizationCases];
            load(filename,'input');
            tester = TestingPhaseFieldHomogenizer(input);
            [xNew,~,~] = tester.compute();

            load(filename,'xRef');
            err = max(norm(xNew(:)-xRef(:))./(norm(xRef(:)+1)));
            tol      = 5e-4;
            testCase.verifyLessThanOrEqual(err, tol)
            close all
        end
        
    end
  
end