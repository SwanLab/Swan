classdef ShFunc_Volume < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        quadrature
        totalVolume
    end

    properties (Access = private)
        mesh
        filter
        designVariableFunction
    end
    
    methods (Access = public)
        function obj = ShFunc_Volume(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function compute(obj)
            obj.computeFunction();
            obj.computeGradient();
            obj.filterGradient();
        end      
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh                   = cParams.mesh;
            obj.filter                 = cParams.filter;
            obj.designVariableFunction = cParams.x;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createTotalVolume(obj)
            q  = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);
            obj.totalVolume = sum(dV(:));
        end

        function computeFunction(obj)
            iPar.mesh     = obj.mesh;
            iPar.quadType = obj.quadrature.order;
            int           = IntegratorFunction(iPar);
            x             = obj.designVariableFunction;
            volume        = int.compute(x);
            obj.value     = volume/(obj.totalVolume);
        end

        function computeGradient(obj)
            f            = AnalyticalFunction.create(@(x) ones(size(x(1,:,:))),1,obj.mesh);
            s.mesh       = obj.mesh;
            s.type       = 'ShapeFunction';
            s.quadType   = 'LINEAR';
            int          = RHSintegrator.create(s);
            test         = P1Function.create(obj.mesh,1);
            Dxv          = int.compute(f,test);
            s.fValues    = Dxv/obj.totalVolume;
            s.mesh       = obj.mesh;
            obj.gradient = P1Function(s);
        end

        function filterGradient(obj)
            g       = obj.gradient;
            regGrad = obj.filter.compute(g,'LINEAR');
            obj.gradient = regGrad;
        end
    end
end

