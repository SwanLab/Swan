classdef ShFunc_InternalEnergySplit < handle
    
    properties (Access = private)
        mesh
        materialPhaseField
    end
    
    methods (Access = public)
        
        function obj = ShFunc_InternalEnergySplit(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunction(obj,u,phi,quadOrder)
            Fbulk = obj.computeEnergyBulk(u,phi,quadOrder);
            Fshear = obj.computeEnergyShear(u,phi,quadOrder);
            F = Fbulk + Fshear;
        end

        function [JSplitu, JSplitphi] = computeGradient(obj,u,phi,quadOrder)
            JSplitu = obj.computeGradientDisplacementSplit(u,phi,quadOrder);
            JSplitphi = obj.computeGradientDamageSplit(u,phi,quadOrder);
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

        function F = computeEnergyBulk(obj,u,phi,quadOrder)
            k = obj.materialPhaseField.getBulkFun(phi,'Interpolated');
            bulkFun = k.*trace(SymGrad(u)).^2;
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(bulkFun);
        end

        function F = computeEnergyShear(obj,u,phi,quadOrder)
            mu = obj.materialPhaseField.getShearFun(phi,'Interpolated');
            strainDev = Deviatoric(SymGrad(u));
            shearFun = mu.*DDP(Voigt(strainDev),Voigt(strainDev));
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = int.compute(shearFun);
        end
        
        function Ju = computeGradientDisplacement(obj,u,phi,quadOrder)
            C = obj.materialPhaseField.setMaterial(phi,'Interpolated');
            sigma = DDP(C,Voigt(SymGrad(u)));
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            C = obj.materialPhaseField.setMaterial(phi,'Jacobian');
            dEnergyFun = DDP(Voigt(SymGrad(u)),DDP(C,Voigt(SymGrad(u))));
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,test);
        end


        function Ju = computeGradientDisplacementSplit(obj,u,phi,quadOrder)
            C = obj.materialPhaseField.setMaterial(phi,'Interpolated');
            sigma = DDP(C,Voigt(SymGrad(u))).evaluate([0;0]);
            sigmaK = zeros(3,1);
            sigmaK(1:2) = sum(sigma(1:2));
            sigmaK(3) = 0;
            sigmaVol = sigmaK*0.5;
            sigmaDev = sigma - sigmaVol;
            k = obj.materialPhaseField.getBulkFun(phi,'Interpolated');
            mu = obj.materialPhaseField.getShearFun(phi,'Interpolated');
            strain = SymGrad(u);
            sigmaBulk = k.*Voigt(Spherical(strain));
            sigmaShear = mu.*Voigt(Deviatoric(strain));
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);

            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            JbulkU = 2.*RHS.compute(sigmaBulk,test);
            JshearU = 2.*RHS.compute(sigmaShear,test);

            Ju = JbulkU + JshearU;
        end






        function Jphi = computeGradientDamageSplit(obj,u,phi,quadOrder)
            k = obj.materialPhaseField.getBulkFun(phi,'Jacobian');
            mu = obj.materialPhaseField.getShearFun(phi,'Jacobian');
            strain = SymGrad(u);
            strainDev = Deviatoric(strain);
            dBulkFun = k.*trace(SymGrad(u)).^2;
            dShearFun = mu.*DDP(Voigt(strainDev),Voigt(strainDev));
            test = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            JbulkPhi = 0.5*RHS.compute(dBulkFun,test);
            JshearPhi = RHS.compute(dShearFun,test);

            Jphi = JbulkPhi + JshearPhi;


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
            ddEnergyFun = DDP(Voigt(SymGrad(u)),DDP(C,Voigt(SymGrad(u))));
            
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