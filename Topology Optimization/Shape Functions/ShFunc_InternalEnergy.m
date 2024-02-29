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
            C = obj.materialPhaseField.setMaterial(phi,'Interpolated');
            energyFun = DDP(SymGrad(u),C,SymGrad(u));
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);
        end

        function [Ju, Jphi] = computeGradient(obj,u,phi,quadOrder)
            Ju = obj.computeGradientDisplacement(u,phi,quadOrder);
            Jphi = obj.computeGradientDamage(u,phi,quadOrder);
        end

        function [Huu, Hphiphi] = computeHessian(obj,u,phi,quadOrder)
            Huu = obj.computeHessianDisplacement(u,phi,quadOrder);
            Hphiphi = obj.computeHessianDamage(u,phi,quadOrder);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialPhaseField = cParams.materialPhaseField;
        end
        
        function Ju = computeGradientDisplacement(obj,u,phi,quadOrder)
            C = obj.materialPhaseField.setMaterial(phi,'Interpolated');
            sigma = DDP(C,SymGrad(u));
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            C = obj.materialPhaseField.setMaterial(phi,'Jacobian');
            dEnergyFun = DDP(SymGrad(u),C,SymGrad(u));
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,test);
        end

        function Huu = computeHessianDisplacement(obj,u,phi,quadOrder)
            mat = obj.materialPhaseField.setMaterial(phi,'Interpolated'); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = u;
            s.material = mat;
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            Huu = LHS.compute();
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)
            C = obj.materialPhaseField.setMaterial(phi,'Hessian');
            ddEnergyFun = DDP(SymGrad(u),C,SymGrad(u));
            
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