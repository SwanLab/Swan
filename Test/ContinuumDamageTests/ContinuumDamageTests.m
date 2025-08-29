classdef ContinuumDamageTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        singleElementCases = {'Linear','Exponential'}
        complexCases = {'SEN'}
    end

    methods (Test, TestTags = {'CD'})
        function testContinuumDamageSingleElem(testCase,singleElementCases)
            filename = ['testContinuumDamageSingleElem',singleElementCases];
            load(filename,'input');
            tester = TestingContinuumDamage(input);
            outputData = tester.compute();
            xNew = outputData.force;
            load(filename,'xRef');
            err = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
        end

        function testContinuumDamageComplexCases(testCase,complexCases)
            filename = ['testContinuumDamage',complexCases];
            load(filename,'input');
            tester = TestingContinuumDamage(input);
            outputData = tester.compute();
            xNew = outputData.force;
            load(filename,'xRef');
            err = norm(xNew-xRef)/norm(xRef);
            tol      = 1e-6;
            testCase.verifyLessThanOrEqual(err, tol)
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