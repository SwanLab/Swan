classdef RHSintegratorConvDifSUSystem < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        function obj = RHSintegratorConvDifSUSystem(cParams)
            obj.init(cParams)
        end

        function f = compute(obj,~,~,source,test)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = 'QUADRATIC';
            int        = RHSintegrator.create(s);
            f          = int.compute(source,test);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
        end
    end
end