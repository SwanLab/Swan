classdef FileComparatorTest < handle

    properties (Access = private)
        testName  
        computerType
        testResultsName
    end

    methods (Access = public)

        function obj = FileComparatorTest(cParams)
            obj.init(cParams);
            obj.computeProblem(cParams);
        end

        function areEqual = computeError(obj)
            refFile      = [obj.testResultsName,'1','.flavia.res'];
            currentFile  = ['Output/',obj.testName,'/',obj.testName,'1','.flavia.res'];
            c            = FileComparator;
            areDifferent = c.areFilesDifferent(refFile,currentFile);
            areEqual     = ~areDifferent;
            delete(currentFile);
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.testName         = cParams.testName;
            obj.computerType     = cParams.computerType;
            obj.testResultsName  = cParams.testName;
            if isfield(cParams, 'testResultsName')
                obj.testResultsName = cParams.testResultsName;
            end
        end

        function computeProblem(obj, s)
            s.testName = obj.testName;
            testComputer = TestComputer.create(obj.computerType, s);
            testComputer.compute();
        end

    end

end