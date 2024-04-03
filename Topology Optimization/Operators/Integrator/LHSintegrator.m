classdef LHSintegrator < handle

    properties (Access = protected)
        mesh
        fun
        quadrature
        quadratureOrder
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = LHSintegratorFactory(); % "Inicialitzem" la classe LHSintegratorFactory
            obj = f.create(s); % Li passem a la classe el que estem resolent (Stokes)
        end

    end

    methods (Access = public)
        
        function obj = LHSintegrator(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end
 
    end

    methods (Access = protected)

        function init(obj, cParams)
            obj.fun      = cParams.fun;
            obj.mesh     = cParams.mesh;
            obj.setQuadratureOrder(cParams);
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = 'QUADRATIC';%obj.fun.order;
            end
        end
        
        function createQuadrature(obj)
            quad = Quadrature.set(obj.mesh.type);
            quadOrder = obj.fun.orderTextual();
            quad.computeQuadrature(quadOrder);
            obj.quadrature = quad;
        end

    end
    
end