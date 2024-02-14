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
            alphaFun = obj.createDissipationFunction(phi);
            
            s.mesh = obj.mesh;
            s.quadType = quadOrder;
            s.type = 'Function';
            int = Integrator.create(s);
            F = (obj.constant/obj.l0)*int.compute(alphaFun);
        end
        
        function J = computeGradient(obj,phi,quadOrder)
            dAlphaFun =  obj.createFirstDerivativeDissipationFunction(phi);
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            J = (obj.constant/obj.l0)*RHS.compute(dAlphaFun, test);
        end
        
        function H = computeHessian(obj,phi,quadOrder)
            ddAlphaFun =  obj.createSecondDerivativeDissipationFunction(phi);
            
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.test = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.function = ddAlphaFun;
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
        
        function alpha = createDissipationFunction(obj,phi)
            s.mesh = obj.mesh;
            s.handleFunction = obj.dissipationInterpolation.fun;
            s.l2function = phi;
            alpha = CompositionFunction(s);
        end
        
        function dAlpha = createFirstDerivativeDissipationFunction(obj,phi)
            s.mesh = obj.mesh;
            s.handleFunction = obj.dissipationInterpolation.dfun;
            s.l2function = phi;
            dAlpha = CompositionFunction(s);
        end
        
        function ddAlpha = createSecondDerivativeDissipationFunction(obj,phi)
            s.mesh = obj.mesh;
            s.handleFunction = obj.dissipationInterpolation.ddfun;
            s.l2function = phi;
            ddAlpha = CompositionFunction(s);
        end
    end
end