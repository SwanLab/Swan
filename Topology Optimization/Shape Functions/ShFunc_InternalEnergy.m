classdef ShFunc_InternalEnergy < handle
    
    properties (Access = private)
        mesh
        materialPhaseField
    end
    
    methods (Access = public)
        
        function obj = ShFunc_InternalEnergy(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunction(obj,u,phi,quadOrder)
            obj.materialPhaseField.setDesignVariable(u,phi,'');
            C = obj.materialPhaseField.obtainTensor();
            energyFun = DDP(Voigt(SymGrad(u)),DDP(C,Voigt(SymGrad(u))));
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);
        end

        function [Ju, Jphi] = computeGradient(obj,u,phi,quadOrder)
            obj.materialPhaseField.setDesignVariable(u,phi,'');                        
            Ju = obj.computeGradientDisplacement(u,quadOrder);
            Jphi = obj.computeGradientDamage(u,phi,quadOrder);
        end

        function [Huu, Hphiphi] = computeHessian(obj,u,phi,quadOrder)
            obj.materialPhaseField.setDesignVariable(u,phi,'');                        
            Huu = obj.computeHessianDisplacement(u,quadOrder);
            Hphiphi = obj.computeHessianDamage(u,phi,quadOrder);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialPhaseField = cParams.materialPhaseField;
        end
        
        function Ju = computeGradientDisplacement(obj,u,quadOrder)
            C = obj.materialPhaseField.obtainTensor();
            sigma = DDP(C,Voigt(SymGrad(u)));
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            dC = obj.materialPhaseField.obtainTensorDerivative();
            dEnergyFun = DDP(Voigt(SymGrad(u)),DDP(dC,Voigt(SymGrad(u))));
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,test);
        end

        function Huu = computeHessianDisplacement(obj,u,quadOrder)
            C = obj.materialPhaseField.obtainTensor(); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = u;
            s.material = C;
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            Huu = LHS.compute();
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)
            ddC = obj.materialPhaseField.obtainTensorSecondDerivative();
            ddEnergyFun = DDP(Voigt(SymGrad(u)),DDP(ddC,Voigt(SymGrad(u))));
            
            s.function = ddEnergyFun;
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.test  = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            Hphiphi = 0.5*LHS.compute();
        end

    end
    
end