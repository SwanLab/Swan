classdef ShFunc_InternalEnergy < handle
    
    properties (Access = private)
        mesh
        materialPhaseField
    end
    
    methods (Access = public)
        
        function obj = ShFunc_InternalEnergy(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunction(obj,u,phi,quad)
            C = obj.materialPhaseField.setMaterial(phi,'Intepolated');
            energyFun = DDP(Voigt(SymGrad(u)),DDP(C,Voigt(SymGrad(u))));
            int = Integrator.create('Function',obj.mesh,quad.order);
            F = 0.5*int.compute(energyFun);
        end

        function [Ju, Jphi] = computeGradient(obj,u,phi,quad)
            Ju = obj.computeGradientDisplacement(u,phi,quad);
            Jphi = obj.computeGradientDamage(u,phi,quad);
        end

        function [Huu, Hphiphi] = computeHessian(obj,u,phi,quad)
            Huu = obj.computeHessianDisplacement(u,phi,quad);
            Hphiphi = obj.computeHessianDamage(u,phi,quad);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialPhaseField = cParams.materialPhaseField;
        end
        
        function Ju = computeGradientDisplacement(obj,u,phi,quad)
            C = obj.materialPhaseField.setMaterial(phi,'Interpolated');
            sigma = DDP(C,Voigt(SymGrad(u)));
            test = LagrangianFunction.create(obj.mesh, 2, 'P1');

            s.mesh = obj.mesh;
            s.quadratureOrder = 'LINEAR';
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quad)
            C = obj.materialPhaseField.setMaterial(phi,'Jacobian');
            dEnergyFun = DDP(Voigt(SymGrad(u)),DDP(C,Voigt(SymGrad(u))));
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quad.order;
            RHS = RHSintegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,test);
        end

        function Huu = computeHessianDisplacement(obj,u,phi,quad)
            mat = obj.materialPhaseField.setMaterial(phi,'Hessian'); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.fun      = u;
            s.material = mat;
            s.quadratureOrder = quad.order;
            LHS = LHSintegrator.create(s);
            Huu = LHS.compute();
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quad)
            C = obj.materialPhaseField.setMaterial(phi,'Hessian');
            ddEnergyFun =  DDP(Voigt(SymGrad(u)),DDP(C,Voigt(SymGrad(u))));
            
            s.function = ddEnergyFun;
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.test  = LagrangianFunction.create(obj.mesh, phi.ndimf, 'P1');
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quad.order;
            LHS = LHSintegrator.create(s);
            Hphiphi = 0.5*LHS.compute();

        end

    end
    
end