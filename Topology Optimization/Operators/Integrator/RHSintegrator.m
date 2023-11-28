classdef RHSintegrator < handle

    properties (Access = protected)
        mesh
        quadrature
        quadratureOrder
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = RHSintegratorFactory();
            obj = f.create(s);
        end

    end
    
    methods (Access = public)
        
%         function obj = RHSintegrator(cParams)
%         end

        function createQuadrature(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature(obj.quadratureOrder);
            obj.quadrature = q;
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = 'LINEAR';
            end
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
        end
        
    end
    
end