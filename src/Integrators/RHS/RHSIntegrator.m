classdef RHSIntegrator < handle

    properties (Access = protected)
        mesh
        quadrature
        quadratureOrder
        test
    end

    methods (Access = public)
        
        function obj = RHSIntegrator(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function RHS = compute(obj,f)
            rhs = obj.computeElementalRHS(f);
            RHS = assembleVector(rhs, obj.test);
        end
    end
    

    methods (Access = public, Static)
        
        function obj = create(s)
            f = RHSIntegratorFactory();
            obj = f.create(s);
        end

    end
    
    methods (Access = public)

        function rhs = computeElementalRHS(obj,f)
            rhs    = zeros(obj.test.nDofsElem,obj.mesh.nelem);
            J      = Jacobian(obj.mesh);
            detJ   = Det(J);
            xV = obj.quadrature.posgp;
            w  = obj.quadrature.weigp;
            v = @(i) Test(obj.test,i);
            for i = 1:obj.test.nDofsElem
                int = (f(v(i)).*detJ)*w';
                rhs(i,:) = rhs(i,:) + squeezeParticular(int.evaluate(xV),2);
            end
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.test = cParams.test;
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end

        function createQuadrature(obj)
            q = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = q;
        end
        
        
    end
    
end