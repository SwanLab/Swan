classdef PhaseFieldTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        singleElementCases = {'Analytic','Split','Homogenized','AnalyticForce'}
        complexCases = {'SEN'}
        homogenizationCases = {'Square','Hexagon','Ellipse'}
    end

    methods (Test, TestTags = {'PF'})
        function testPhaseFieldSingleElem(testCase,singleElementCases)
            filename = ['testPhaseFieldSingleElem',singleElementCases];
            load(filename,'input');
            tester = TestingPhaseField(input);
            outputData = tester.compute();
            xNew = outputData.reaction;
            load(filename,'xRef');
            err = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testPhaseFieldComplexCases(testCase,complexCases)
            filename = ['testPhaseField',complexCases];
            load(filename,'input');
            tester = TestingPhaseField(input);
            outputData = tester.compute();
            xNew = outputData.reaction;
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
            err = max(pagenorm(xNew-xRef)./(pagenorm(xRef)+1));
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
            close all
        end
        
    end

    methods(Static,Access = private)

        function computeTotalEnergy(energyCell)
            E = struct2cell(energyCell);
            totE = zeros(size(E{1}));
            for i=1:length(E)
                totE = totE + E{i};
            end
        end

    end
  
end