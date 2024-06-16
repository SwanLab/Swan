classdef RHSintegratorConvDifSUPGSystem < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        function obj = RHSintegratorConvDifSUPGSystem(cParams)
            obj.init(cParams)
        end

        function f = compute(obj,a,nu,source,test)
            t  = obj.computeStabParameter(a,nu);
            f1 = obj.computeRHSShape(source,test);
            f2 = obj.computeRHSShapeDer(source,test);
            f  = f1 + t*a*f2;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
        end

        function tau = computeStabParameter(obj,a,nu)
            h    = obj.mesh.computeMeanCellSize();
            Pe   = a*h/(2*nu);
            alfa = coth(Pe)-1/Pe;
            tau  = alfa*h/(2*a);
        end

        function f = computeRHSShape(obj,source,test)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = 'QUADRATIC';
            int        = RHSintegrator.create(s);
            f          = int.compute(source,test);
        end

        function f = computeRHSShapeDer(obj,source,test)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadType = 'QUADRATIC';
            int        = RHSintegrator.create(s);
            f          = int.compute(source,test);
        end

    end
end