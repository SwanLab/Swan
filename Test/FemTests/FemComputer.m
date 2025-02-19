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

    methods (Access = public) % Aquí arribem desde "testComputer = TestComputer.create(obj.computerType, s);", a PrecomputedVariableTest, passant per TestComputer.
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
            d_u = obj.computation.uFun.fValues;
            obj.variables.d_u = reshape(d_u', [numel(d_u) 1]);
        end
    end

end