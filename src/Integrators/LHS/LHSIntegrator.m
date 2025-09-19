classdef LHSIntegrator < handle

    properties (Access = protected)
        mesh
        test, trial
        quadratureOrder        
        quadrature
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = LHSIntegratorFactory();
            obj = f.create(s);
        end

    end

    methods (Access = public)
        
        function obj = LHSIntegrator(cParams)
            obj.init(cParams);
            obj.createQuadrature();
        end

        function LHS = compute(obj,f)
            lhs = obj.computeElementalLHS(f);
            LHS = assembleMatrix(lhs, obj.test, obj.trial);
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj,f)
            lhs    = zeros(obj.test.nDofsElem,obj.trial.nDofsElem,obj.mesh.nelem);
            detJ   = DetJacobian(obj.mesh);
            %detJ   = Det(J);
            xV = obj.quadrature.posgp;
            w  = obj.quadrature.weigp;
            v = @(i) Test(obj.test,i);
            u = @(j) Test(obj.trial,j);              
            for i = 1:obj.test.nDofsElem
                for j = 1:obj.trial.nDofsElem
                    int = (f(u(j),v(i)).*detJ)*w';
                    lhs(i,j,:) = lhs(i,j,:) + int.evaluate(xV);
                end
            end
        end

        function init(obj, cParams)
            obj.test  = cParams.test;
            obj.trial = cParams.trial;
            obj.mesh  = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end
    
        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh, obj.quadratureOrder);
            obj.quadrature = quad;
        end

    end
    
end