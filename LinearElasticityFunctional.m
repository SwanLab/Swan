classdef LinearElasticityFunctional < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        material
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LinearElasticityFunctional (cParams)
            obj.init(cParams)
        end

        function val = compute(obj, uFun)
            val = 1;
        end
        

        function Ju = computeGradient(obj, uFun)
            strainVgt = SymGrad(uFun);
%             strainVgt.ndimf = 4;
            sigma = DDP(obj.material,strainVgt);
%             sigma.ndimf = 3;
            test = LagrangianFunction.create(obj.mesh, uFun.ndimf, uFun.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = 3;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Huu = computeHessian(obj, uFun) 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.material = obj.material;
            s.quadratureOrder = 3;
            s.test     = LagrangianFunction.create(obj.mesh,uFun.ndimf, 'P1');
            s.trial    = uFun;
            LHS = LHSintegrator.create(s);
            Huu = LHS.compute();
        end


    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.material = cParams.material;
        end

        
    end
    
end