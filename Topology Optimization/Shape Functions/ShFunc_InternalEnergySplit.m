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

        function dEu = computeGradientDisplacement(obj,u,phi,quadOrder)
            dEvol = obj.computeVolumetricEnergyDisplacementGradient(u,phi,quadOrder);
            dEdev = obj.computeDeviatoricEnergyDisplacementGradient(u,phi,quadOrder);
            dEu = dEvol + dEdev;
        end

        function dEphi = computeGradientDamage(obj,u,phi,quadOrder)
            dEvol = obj.computeVolumetricEnergyDamageGradient(u,phi,quadOrder);
            dEdev = obj.computeDeviatoricEnergyDamageGradient(u,phi,quadOrder);
            dEphi = dEvol + dEdev;
        end

        function ddEuu = computeHessianDisplacement(obj,u,phi,quadOrder)
           ddEvol = obj.computeVolumetricEnergyDisplacementHessian(u,phi,quadOrder);
           ddEdev = obj.computeDeviatoricEnergyDisplacementHessian(u,phi,quadOrder);
           ddEuu = ddEvol + ddEdev;
        end

        function ddEphiphi = computeHessianDamage(obj,u,phi,quadOrder)        
            ddEvol  = obj.computeVolumetricEnergyDamageHessian(u,phi,quadOrder);
            ddEdev   = obj.computeDeviatoricEnergyDamageHessian(u,phi,quadOrder);
            ddEphiphi = ddEvol + ddEdev;
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialPhaseField = cParams.material;
        end

        function F = computeEnergyBulk(obj,u,phi,quadOrder)
            k = obj.materialPhaseField.getBulkFun(u,phi,'Interpolated');
            bulkFun = VolumetricElasticEnergyDensity(u,k);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = int.compute(bulkFun);            
        end

        function F = computeEnergyShear(obj,u,phi,quadOrder)
            mu = obj.materialPhaseField.getShearFun(u,phi,'Interpolated');
            shearFun = DeviatoricElasticEnergyDensity(u,mu);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = int.compute(shearFun);            
        end

        
        function dE = computeVolumetricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            Cbulk    = obj.materialPhaseField.getBulkMaterial(u,phi,'Interpolated');
            sigmaVol = DDP(Cbulk,SymGrad(u));
            dE  = obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaVol,u,quadOrder);
        end

        function dE = computeDeviatoricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            Cshear = obj.materialPhaseField.getShearMaterial(u,phi,'Interpolated');
            sigmaDev = DDP(Cshear,SymGrad(u));
            dE =  obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaDev,u,quadOrder);
        end        

        function F = computeShapeSymmetricDerivativeIntegralWithField(obj,f,u,quadOrder)
            test = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSintegrator.create(s);
            F = RHS.compute(f,test);            
        end 

        function dE = computeVolumetricEnergyDamageGradient(obj,u,phi,quadOrder)
            dk    = obj.materialPhaseField.getBulkFun(u,phi,'Jacobian');
            deVol = VolumetricElasticEnergyDensity(u,dk);
            dE    =  obj.computeShapeIntegralWithField(deVol,phi,quadOrder);
        end

        function dE = computeDeviatoricEnergyDamageGradient(obj,u,phi,quadOrder)
            dmu = obj.materialPhaseField.getShearFun(u,phi,'Jacobian');
            deDev = DeviatoricElasticEnergyDensity(u,dmu);
            dE =  obj.computeShapeIntegralWithField(deDev,phi,quadOrder);
        end        

        function F = computeShapeIntegralWithField(obj,f,phi,quadOrder)
            test   = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSintegrator.create(s);
            F   =  RHS.compute(f,test);
        end        


        function ddE = computeVolumetricEnergyDisplacementHessian(obj,u,phi,quadOrder)
            Cbulk = obj.materialPhaseField.getBulkMaterial(u,phi,'Interpolated');
            ddE   = obj.computeElasticStiffnesMatrix(Cbulk,u,quadOrder);
        end

        function ddE = computeDeviatoricEnergyDisplacementHessian(obj,u,phi,quadOrder)
            Cshear = obj.materialPhaseField.getShearMaterial(u,phi,'Interpolated');
            ddE    = obj.computeElasticStiffnesMatrix(Cshear,u,quadOrder);
        end        

        function K = computeElasticStiffnesMatrix(obj,C,u,quadOrder)
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.test  = LagrangianFunction.create(obj.mesh, u.ndimf, u.order);
            s.material = C;
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            K = LHS.compute();
        end        

        function ddE = computeVolumetricEnergyDamageHessian(obj,u,phi,quadOrder)
            ddk    = obj.materialPhaseField.getBulkFun(u,phi,'Hessian');
            ddeVol = VolumetricElasticEnergyDensity(u,ddk);
            ddE    = obj.computeMassWithFunction(ddeVol,phi,quadOrder);
        end

        function ddE = computeDeviatoricEnergyDamageHessian(obj,u,phi,quadOrder)
            ddmu   = obj.materialPhaseField.getShearFun(u,phi,'Hessian');
            ddeDev = DeviatoricElasticEnergyDensity(u,ddmu);
            ddE    = obj.computeMassWithFunction(ddeDev,phi,quadOrder);
        end

        function Mf = computeMassWithFunction(obj,f,phi,quadOrder)
            s.fun = f;
            s.trial = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.test  = LagrangianFunction.create(obj.mesh, phi.ndimf, phi.order);
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSintegrator.create(s);
            Mf = LHS.compute();
        end        
        
    end
    
end