classdef FemComputer < handle

    properties (Access = public)
        computation
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
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
            s.interpolationType = obj.interpolationType;
            obj.computation = FEM.create(s);
            obj.computation.solve();
        end
    end

end