classdef PhaseFieldInternalEnergySplitFunctional < handle
    
    properties (Access = private)
        mesh
        materialPhaseField
        testU
        testPhi
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldInternalEnergySplitFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeCost(obj,u,phi,quadOrder)
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
            mPF     = obj.materialPhaseField;
            k       = mPF.obtainBulkFunction(u,phi);
            bulkFun = VolumetricElasticEnergyDensity(u,k);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = int.compute(bulkFun);            
        end

        function F = computeEnergyShear(obj,u,phi,quadOrder)
            mPF      = obj.materialPhaseField;
            mu       = mPF.obtainShearFunction(phi);
            shearFun = DeviatoricElasticEnergyDensity(u,mu);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = int.compute(shearFun);            
        end


        function dE = computeVolumetricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            mPF      = obj.materialPhaseField;
            Cbulk    = mPF.obtainTensorVolumetric(u,phi);
            sigmaVol = DDP(Cbulk,SymGrad(u));
            dE  = obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaVol,quadOrder);
        end

        function dE = computeDeviatoricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            mPF      = obj.materialPhaseField;
            Cshear   = mPF.obtainTensorDeviatoric(phi);
            sigmaDev = DDP(Cshear,SymGrad(u));
            dE =  obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaDev,quadOrder);
        end        

        function F = computeShapeSymmetricDerivativeIntegralWithField(obj,f,quadOrder)
            F = IntegrateRHS(@(v) DDP(SymGrad(v),f),obj.testU,obj.mesh,'Domain',quadOrder);
        end 


        function dE = computeVolumetricEnergyDamageGradient(obj,u,phi,quadOrder)
            mPF   = obj.materialPhaseField;
            dk    = mPF.obtainBulkDerivative(u,phi);
            deVol = VolumetricElasticEnergyDensity(u,dk);
            dE    = IntegrateRHS(@(v) DP(v,deVol),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end

        function dE = computeDeviatoricEnergyDamageGradient(obj,u,phi,quadOrder)
            mPF   = obj.materialPhaseField;
            dmu   = mPF.obtainShearDerivative(phi);
            deDev = DeviatoricElasticEnergyDensity(u,dmu);
            dE    = IntegrateRHS(@(v) DP(v,deDev),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end        

        function ddE = computeVolumetricEnergyDisplacementHessian(obj,u,phi,quadOrder)
            mPF   = obj.materialPhaseField;
            Cbulk = mPF.obtainTensorVolumetric(u,phi);
            ddE   = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(Cbulk,SymGrad(u))),obj.testU,obj.testU,obj.mesh,'Domain',quadOrder);
        end

        function ddE = computeDeviatoricEnergyDisplacementHessian(obj,phi,quadOrder)
            mPF    = obj.materialPhaseField;
            Cshear = mPF.obtainTensorDeviatoric(phi);
            ddE    = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(Cshear,SymGrad(u))),obj.testU,obj.testU,obj.mesh,'Domain',quadOrder);
        end        

        function ddE = computeVolumetricEnergyDamageHessian(obj,u,phi,quadOrder)
            mPF    = obj.materialPhaseField;
            ddk    = mPF.obtainBulkSecondDerivative(u,phi);
            ddeVol = VolumetricElasticEnergyDensity(u,ddk);
            ddE    = obj.computeMassWithFunction(ddeVol,quadOrder);
        end

        function ddE = computeDeviatoricEnergyDamageHessian(obj,u,phi,quadOrder)
            mPF    = obj.materialPhaseField;
            ddmu   = mPF.obtainShearSecondDerivative(phi);
            ddeDev = DeviatoricElasticEnergyDensity(u,ddmu);
            ddE    = obj.computeMassWithFunction(ddeDev,quadOrder);
        end

        function Mf = computeMassWithFunction(obj,f,quadOrder)
            Mf = IntegrateLHS(@(u,v) f.*DP(v,u),obj.testPhi,obj.testPhi,obj.mesh,'Domain',quadOrder);
        end        
        
    end
    
end