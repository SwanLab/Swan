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
            refFile      = [obj.testResultsName,'_1','.vtu'];
            currentFile  = ['Output/',obj.testName,'/',obj.testName,'_1','.vtu'];
            c            = FileComparator;
            areDifferent = c.areFilesDifferent(currentFile,refFile);
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