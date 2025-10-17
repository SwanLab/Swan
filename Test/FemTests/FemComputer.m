classdef FemComputer < handle

    properties (Access = public)
        computation
        variables
        solverType
    end

    properties (Access = private)
        testName
        interpolationType
    end

    methods (Access = public)
        function obj = FemComputer(cParams)
            obj.testName          = cParams.testName;
            if isfield(cParams, 'interpolationType')
                obj.interpolationType = cParams.interpolationType;
            else
                obj.interpolationType = 'LINEAR';
            end

            if isfield(cParams, 'solverType')
                obj.solverType = cParams.solverType;
            else
                obj.solverType = 'REDUCED';
            end
            
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
            s.interpolationType = obj.interpolationType;
            s.solverType = obj.solverType;
            s.solverMode = 'DISP';
            obj.computation = PhysicalProblem.create(s);
            obj.computation.solve();
            obj.variables.d_u = obj.computation.uFun.fValues;
        end
    end

end