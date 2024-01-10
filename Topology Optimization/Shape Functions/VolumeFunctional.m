classdef VolumeFunctional < handle
    
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
    end
    
    methods (Access = public)
        function obj = VolumeFunctional(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createTotalVolume();
        end

        function computeFunctionAndGradient(obj,x)
            obj.computeFunction(x);
            obj.computeGradient(x);
        end      
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('LINEAR');
            obj.quadrature = quad;
        end

        function createTotalVolume(obj)
            dV = obj.mesh.computeDvolume(obj.quadrature);
            obj.totalVolume = sum(dV(:));
        end

        function computeFunction(obj,x)
            s.mesh     = obj.mesh;
            s.quadType = obj.quadrature.order;
            int        = IntegratorFunction(s);
            volume     = int.compute(x);
            obj.value  = volume/obj.totalVolume;
        end

        function computeGradient(obj,x)
            s.mesh      = obj.mesh;
            s.feFunType = class(x);
            s.ndimf     = 1;
            g           = FeFunction.createEmpty(s); % create
            g.fValues   = ones(x.nDofs,1)/obj.totalVolume;
            obj.gradient = g;
        end
    end
end

