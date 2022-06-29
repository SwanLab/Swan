classdef storedVariableTest < handle

    properties (Access = protected)
        testName  
        computerType
        testResultsName
        nIter
    end

    methods (Access = public)

        function obj = storedVariableTest(cParams)
            obj.init(cParams);
            obj.computeProblem(cParams);
        end

        function areEqual = computeError(obj)
            refFile     = [obj.testResultsName,num2str(obj.nIter),'.flavia.res'];
            currentFile = ['Output/',obj.testName,'/',obj.testName,num2str(obj.nIter),'.flavia.res'];
            compare     = FileComparator;
            areEqual    = ~compare.areFilesDifferent(refFile,currentFile);
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
            obj.nIter =  testComputer.computation.optimizer.nIter - 1;
        end

    end

end