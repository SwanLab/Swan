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
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            J = (obj.constant/obj.l0)*RHS.compute(dAlphaFun, test);
        end
        
        function H = computeHessian(obj,phi,quadOrder)
            ddAlphaFun =  obj.createDissipationFunction(phi,'Hessian');
            
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.test = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.fun = ddAlphaFun;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            H = (obj.constant/obj.l0)*LHS.compute();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.dissipationInterpolation = cParams.dissipationInterpolation;
            obj.constant = cParams.constant;
            obj.l0 = cParams.l0;
        end

        function alpha = createDissipationFunction(obj,phi,type)
            switch type
                case 'Function'
                    s.handleFunction = obj.dissipationInterpolation.fun;
                case'Jacobian'
                    s.handleFunction = obj.dissipationInterpolation.dfun;
                case 'Hessian'
                    s.handleFunction = obj.dissipationInterpolation.ddfun;
            end
            s.mesh = obj.mesh;
            s.l2function = phi;
            alpha = CompositionFunction(s);
        end

    end
end