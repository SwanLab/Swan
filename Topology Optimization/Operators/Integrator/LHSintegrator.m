classdef LHSintegrator < handle

    properties (Access = protected)
        mesh
        fun
        quadrature
        quadratureOrder
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = LHSintegratorFactory();
            obj = f.create(s);
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
                warning('Assuming quadrature order = 2')
                obj.quadratureOrder = 2;%obj.fun.order;
            end
        end
        
        function createQuadrature(obj)
%             quadOrder = obj.fun.getOrderNum();
            quad = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = quad;
        end

    end
    
end