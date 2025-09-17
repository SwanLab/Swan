classdef RHSIntegrator < handle

    properties (Access = protected)
        mesh
        quadrature
        quadratureOrder
        test
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = RHSIntegratorFactory();
            obj = f.create(s);
        end

        function obj = compute(obj,f)

        end

    end
    
    methods (Access = public)

        function createQuadrature(obj)
            q = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = q;
        end

        function setQuadratureOrder(obj, cParams)
            if isfield(cParams, 'quadratureOrder')
                obj.quadratureOrder = cParams.quadratureOrder;
            else
                obj.quadratureOrder = 2;
            end
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.test = cParams.test;
            obj.mesh = cParams.mesh;
            obj.setQuadratureOrder(cParams);            
        end
        
    end
    
end