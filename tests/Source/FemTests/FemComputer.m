classdef FemComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = FemComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
            obj.computation = FEM.create(s);
            obj.computation.solve();
        end
    end

end