classdef ShFunc_LocalDamage < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        mesh
        dissipationInterpolation
        constant
        l0
        test
    end
    
    methods (Access = public)
        
        function obj = ShFunc_LocalDamage(cParams)
            obj.init(cParams)
        end
        
        function F = computeFunction(obj,phi,quadOrder)
            alphaFun = obj.createDissipationFunction(phi,'Function');
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = (obj.constant/obj.l0)*int.compute(alphaFun);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            dAlphaFun =  obj.createDissipationFunction(phi,'Jacobian');
            test = obj.test;
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSIntegrator.create(s);
            J = (obj.constant/obj.l0)*RHS.compute(dAlphaFun, test);
        end
        
        function H = computeHessian(obj,phi,quadOrder)
            ddAlphaFun =  obj.createDissipationFunction(phi,'Hessian');
            
            s.trial = obj.test;
            s.test = obj.test;
            s.fun = ddAlphaFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSIntegrator.create(s);
            H = (obj.constant/obj.l0)*LHS.compute();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.test = LagrangianFunction.create(obj.mesh, 1, 'P1');            
            obj.dissipationInterpolation = cParams.dissipation.interpolation;
            obj.constant = cParams.dissipation.constant;
            obj.l0 = cParams.l0;
        end

        function disFun = createDissipationFunction(obj,phi,type)
            switch type
                case 'Function'
                    fun = obj.dissipationInterpolation.fun;
                case'Jacobian'
                     fun = obj.dissipationInterpolation.dfun;
                case 'Hessian'
                    fun = obj.dissipationInterpolation.ddfun;
            end
            s.operation = @(xV) fun.evaluate(phi.evaluate(xV));
            s.ndimf = 1;
            s.mesh  = obj.mesh;
            disFun = DomainFunction(s);
        end

    end
end