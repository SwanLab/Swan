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
            obj.materialPhaseField = cParams.material;
        end

        function F = computeEnergyBulk(obj,u,phi,quadOrder)
            k = obj.materialPhaseField.getBulkFun(u,phi,'Interpolated');
            bulkFun = k.*trace(SymGrad(u)).^2;
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(bulkFun);
        end

        function F = computeEnergyShear(obj,u,phi,quadOrder)
            mu = obj.materialPhaseField.getShearFun(u,phi,'Interpolated');
            strainDev = Deviatoric(SymGrad(u));
            shearFun = mu.*DDP(Voigt(strainDev),Voigt(strainDev));
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = int.compute(shearFun);
        end
        
        function Ju = computeGradientDisplacement(obj,u,phi,quadOrder)
            %k = obj.materialPhaseField.getBulkFun(u,phi,'Interpolated');
            Cbulk = obj.materialPhaseField.getBulkMaterial(u,phi);
            %mu = obj.materialPhaseField.getShearFun(u,phi,'Interpolated');
            Cshear = obj.materialPhaseField.getShearMaterial(u,phi);
            strain = SymGrad(u);
            %sigmaBulk = 2.*k.*VoigtStress(Spherical(strain)); % 2 is from ndim = 2;
            %sigmaShear = 2.*mu.*VoigtStress(Deviatoric(strain));
            sigmaBulk = DDP(Cbulk,Voigt(strain));
            sigmaShear = DDP(Cshear,Voigt(strain));


            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            JbulkU = RHS.compute(sigmaBulk,test);
            JshearU = RHS.compute(sigmaShear,test);

            Ju = JbulkU + JshearU;
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            k = obj.materialPhaseField.getBulkFun(u,phi,'Jacobian');
            mu = obj.materialPhaseField.getShearFun(u,phi,'Jacobian');
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
            Cbulk = obj.materialPhaseField.getBulkMaterial(u,phi);
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.test  = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.material = Cbulk;
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            HbulkUU = LHS.compute();

            Cshear = obj.materialPhaseField.getShearMaterial(u,phi);
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.test  = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.material = Cshear;
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            HshearUU = LHS.compute();

            Huu = HbulkUU + HshearUU;
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)
            strain = SymGrad(u);
            strainDev = Deviatoric(strain);

            k = obj.materialPhaseField.getBulkFun(u,phi,'Hessian');
            ddBulkFun = k.*trace(SymGrad(u)).^2;
            s.fun = ddBulkFun;
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.test  = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            HbulkPhiPhi = 0.5*LHS.compute();

            mu = obj.materialPhaseField.getShearFun(u,phi,'Hessian');
            ddShearFun = mu.*DDP(Voigt(strainDev),Voigt(strainDev));
            s.function = ddShearFun;
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.test  = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            HshearPhiPhi = LHS.compute();

            Hphiphi = HbulkPhiPhi + HshearPhiPhi;
        end
        
    end
    
end