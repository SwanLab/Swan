classdef InternalEnergySplitFunctional < handle
    
    properties (Access = private)
        mesh
        materialPhaseField
        testU
        testPhi
    end
    
    methods (Access = public)
        
        function obj = InternalEnergySplitFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunctional(obj,u,phi,quadOrder)
            Fbulk  = obj.computeEnergyBulk(u,phi,quadOrder);
            Fshear = obj.computeEnergyShear(u,phi,quadOrder);
            F      = Fbulk + Fshear;
        end

        function dEu = computeGradientDisplacement(obj,u,phi,quadOrder)
            dEvol = obj.computeVolumetricEnergyDisplacementGradient(u,phi,quadOrder);
            dEdev = obj.computeDeviatoricEnergyDisplacementGradient(u,phi,quadOrder);
            dEu   = dEvol + dEdev;
        end

        function dEphi = computeGradientDamage(obj,u,phi,quadOrder)
            dEvol = obj.computeVolumetricEnergyDamageGradient(u,phi,quadOrder);
            dEdev = obj.computeDeviatoricEnergyDamageGradient(u,phi,quadOrder);
            dEphi = dEvol + dEdev;
        end

        function ddEuu = computeHessianDisplacement(obj,u,phi,quadOrder)
           ddEvol = obj.computeVolumetricEnergyDisplacementHessian(u,phi,quadOrder);
           ddEdev = obj.computeDeviatoricEnergyDisplacementHessian(phi,quadOrder);
           ddEuu  = ddEvol + ddEdev;
        end

        function ddEphiphi = computeHessianDamage(obj,u,phi,quadOrder)        
            ddEvol    = obj.computeVolumetricEnergyDamageHessian(u,phi,quadOrder);
            ddEdev    = obj.computeDeviatoricEnergyDamageHessian(u,phi,quadOrder);
            ddEphiphi = ddEvol + ddEdev;
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.materialPhaseField = cParams.material;
            obj.testPhi  = copy(cParams.testSpace.phi);
            obj.testU    = copy(cParams.testSpace.u);
        end

        function F = computeEnergyBulk(obj,u,phi,quadOrder)
            k       = obj.materialPhaseField.getBulkFun(u,phi);
            bulkFun = VolumetricElasticEnergyDensity(u,k);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = int.compute(bulkFun);            
        end

        function F = computeEnergyShear(obj,u,phi,quadOrder)
            mu       = obj.materialPhaseField.getShearFun(phi,'Interpolated');
            shearFun = DeviatoricElasticEnergyDensity(u,mu);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = int.compute(shearFun);            
        end


        function dE = computeVolumetricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            Cbulk    = obj.materialPhaseField.getBulkMaterial(u,phi,'Interpolated');
            sigmaVol = DDP(Cbulk,SymGrad(u));
            dE  = obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaVol,quadOrder);
        end

        function dE = computeDeviatoricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            Cshear   = obj.materialPhaseField.getShearMaterial(phi,'Interpolated');
            sigmaDev = DDP(Cshear,SymGrad(u));
            dE =  obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaDev,quadOrder);
        end        

        function F = computeShapeSymmetricDerivativeIntegralWithField(obj,f,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder = quadOrder;
            RHS = RHSIntegrator.create(s);
            F   = RHS.compute(f,obj.testU);            
        end 


        function dE = computeVolumetricEnergyDamageGradient(obj,u,phi,quadOrder)
            dk    = obj.materialPhaseField.getBulkFun(u,phi,'Jacobian');
            deVol = VolumetricElasticEnergyDensity(u,dk);
            dE    =  obj.computeShapeIntegralWithField(deVol,quadOrder);
        end

        function dE = computeDeviatoricEnergyDamageGradient(obj,u,phi,quadOrder)
            dmu   = obj.materialPhaseField.getShearFun(phi,'Jacobian');
            deDev = DeviatoricElasticEnergyDensity(u,dmu);
            dE    =  obj.computeShapeIntegralWithField(deDev,quadOrder);
        end        

        function F = computeShapeIntegralWithField(obj,f,quadOrder)        
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSIntegrator.create(s);
            F   =  RHS.compute(f,obj.testPhi);
        end        


        function ddE = computeVolumetricEnergyDisplacementHessian(obj,u,phi,quadOrder)
            Cbulk = obj.materialPhaseField.getBulkMaterial(u,phi,'Interpolated');
            ddE   = obj.computeElasticStiffnesMatrix(Cbulk,quadOrder);
        end

        function ddE = computeDeviatoricEnergyDisplacementHessian(obj,phi,quadOrder)
            Cshear = obj.materialPhaseField.getShearMaterial(phi,'Interpolated');
            ddE    = obj.computeElasticStiffnesMatrix(Cshear,quadOrder);
        end        

        function K = computeElasticStiffnesMatrix(obj,C,quadOrder)
            s.trial    = obj.testU;
            s.test     = obj.testU;
            s.material = C;
            s.mesh     = obj.mesh;
            s.type     = 'ElasticStiffnessMatrix';
            s.quadratureOrder = quadOrder;
            LHS = LHSIntegrator.create(s);
            K   = LHS.compute();
        end        

        function ddE = computeVolumetricEnergyDamageHessian(obj,u,phi,quadOrder)
            ddk    = obj.materialPhaseField.getBulkFun(u,phi,'Hessian');
            ddeVol = VolumetricElasticEnergyDensity(u,ddk);
            ddE    = obj.computeMassWithFunction(ddeVol,quadOrder);
        end

        function ddE = computeDeviatoricEnergyDamageHessian(obj,u,phi,quadOrder)
            ddmu   = obj.materialPhaseField.getShearFun(phi,'Hessian');
            ddeDev = DeviatoricElasticEnergyDensity(u,ddmu);
            ddE    = obj.computeMassWithFunction(ddeDev,quadOrder);
        end

        function Mf = computeMassWithFunction(obj,f,quadOrder)
            s.fun   = f;
            s.trial = obj.testPhi;
            s.test  = obj.testPhi;
            s.mesh  = obj.mesh;
            s.type  = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSIntegrator.create(s);
            Mf  = LHS.compute();
        end        
        
    end
    
end