classdef ReadingFilesTests < handle & matlab.unittest.TestCase

    properties (TestParameter)
        readingTests = {'testReadingGmsh'}
    end

    methods (Test, TestTags = {'ReadingFiles', 'Fast', 'Gmsh'})

        function testReading(testCase, readingTests)
            s.computerType    = 'GMSH';
            s.testName         = readingTests;
            s.variablesToStore = {'connec', 'coord', 'isElemInThisSet', ...
                                  'masterSlave', 'corners'};
            test = PrecomputedVariableTest(s);
            err = test.computeError();
            tol = 1e-13;
            testCase.verifyLessThanOrEqual(err, tol)
        end

    end
    
end