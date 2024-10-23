classdef ShFunc_InternalEnergy < handle
    
    properties (Access = private)
        mesh
        material
    end
    
    methods (Access = public)
        
        function obj = ShFunc_InternalEnergy(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunction(obj,u,phi,quadOrder)
            obj.material.setDesignVariable(u,phi);
            C = obj.material.obtainTensor();
            energyFun = DDP(SymGrad(u),DDP(C{1},SymGrad(u)));
           % figure(100); energyFun.project('P1',u.mesh).plot;
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);
        end

        function [Ju, Jphi] = computeGradient(obj,u,phi,quadOrder)
            obj.material.setDesignVariable(u,phi);                        
            Ju = obj.computeGradientDisplacement(u,quadOrder);
            Jphi = obj.computeGradientDamage(u,phi,quadOrder);
        end

        function [Huu, Hphiphi] = computeHessian(obj,u,phi,quadOrder)
            obj.material.setDesignVariable(u,phi);                        
            Huu = obj.computeHessianDisplacement(u,quadOrder);
            Hphiphi = obj.computeHessianDamage(u,phi,quadOrder);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.material = cParams.material;
        end
        
        function Ju = computeGradientDisplacement(obj,u,quadOrder)
            C = obj.material.obtainTensor();
            sigma = DDP(C{1},SymGrad(u));
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            dC = obj.material.obtainTensorDerivative();
            dEnergyFun = DDP(SymGrad(u),DDP(dC{1},SymGrad(u)));
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,test);
        end

        function Huu = computeHessianDisplacement(obj,u,quadOrder)
            C = obj.material.obtainTensor(); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.trial     = u;
            s.test      = u;
            s.material = C{1};
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            Huu = LHS.compute();
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)
            ddC = obj.material.obtainTensorSecondDerivative();
            ddEnergyFun = DDP(SymGrad(u),DDP(ddC{1},SymGrad(u)));
            
            s.fun = ddEnergyFun;
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