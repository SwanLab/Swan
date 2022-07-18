classdef PostProcTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        fastDisp = {
            'test_postprocess'
            }
    end
        
    methods (Test, TestTags = {'PostProc', 'Various', 'FastDisp'})

        function testFastDisplacement(testCase, fastDisp)
            s.computerType    = 'TOPOPT';
            s.testName         = fastDisp;
            s.testResultsName  = [fastDisp,'Ref'];
            test = FileComparatorTest(s);
            err = test.computeError();
            testCase.verifyTrue(err)
        end

    end

end