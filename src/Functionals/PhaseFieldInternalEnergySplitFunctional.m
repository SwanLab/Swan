classdef PhaseFieldInternalEnergySplitFunctional < handle
    
    properties (Access = private)
        mesh
        mu, dmu, d2mu
        k, dk, d2k
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
            obj.mu   = cParams.mu;   obj.dmu  = cParams.dmu;  obj.d2mu = cParams.d2mu;
            obj.k    = cParams.k;    obj.dk   = cParams.dk;   obj.d2k  = cParams.d2k;
            obj.testPhi  = copy(cParams.testSpace.phi);
            obj.testU    = copy(cParams.testSpace.u);
        end

        function F = computeEnergyBulk(obj,u,phi,quadOrder)
            k       = obj.k(phi.fun);
            bulkFun = VolumetricElasticEnergyDensity(u,k);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = int.compute(bulkFun);            
        end

        function F = computeEnergyShear(obj,u,phi,quadOrder)
            mu       = obj.mu(phi.fun);
            shearFun = DeviatoricElasticEnergyDensity(u,mu);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F   = int.compute(shearFun);            
        end


        function dE = computeVolumetricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            Cbulk    = obj.obtainTensorVolumetric(phi);
            sigmaVol = DDP(Cbulk,SymGrad(u));
            dE  = obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaVol,quadOrder);
        end

        function dE = computeDeviatoricEnergyDisplacementGradient(obj,u,phi,quadOrder)
            Cshear   = obj.obtainTensorDeviatoric(phi);
            sigmaDev = DDP(Cshear,SymGrad(u));
            dE =  obj.computeShapeSymmetricDerivativeIntegralWithField(sigmaDev,quadOrder);
        end        

        function F = computeShapeSymmetricDerivativeIntegralWithField(obj,f,quadOrder)
            F = IntegrateRHS(@(v) DDP(SymGrad(v),f),obj.testU,obj.mesh,'Domain',quadOrder);
        end 


        function dE = computeVolumetricEnergyDamageGradient(obj,u,phi,quadOrder)
            dk    = obj.dk(phi.fun);
            deVol = VolumetricElasticEnergyDensity(u,dk);
            dE    = IntegrateRHS(@(v) DP(v,deVol),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end

        function dE = computeDeviatoricEnergyDamageGradient(obj,u,phi,quadOrder)
            dmu   = obj.dmu(phi.fun);
            deDev = DeviatoricElasticEnergyDensity(u,dmu);
            dE    = IntegrateRHS(@(v) DP(v,deDev),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end        

        function ddE = computeVolumetricEnergyDisplacementHessian(obj,u,phi,quadOrder)
            Cbulk = obj.obtainTensorVolumetric(phi);
            ddE   = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(Cbulk,SymGrad(u))),obj.testU,obj.testU,obj.mesh,'Domain',quadOrder);
        end

        function ddE = computeDeviatoricEnergyDisplacementHessian(obj,phi,quadOrder)
            Cshear = obj.obtainTensorDeviatoric(phi);
            ddE    = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(Cshear,SymGrad(u))),obj.testU,obj.testU,obj.mesh,'Domain',quadOrder);
        end        

        function ddE = computeVolumetricEnergyDamageHessian(obj,u,phi,quadOrder)
            ddk    = obj.d2k(phi.fun);
            ddeVol = VolumetricElasticEnergyDensity(u,ddk);
            ddE    = obj.computeMassWithFunction(ddeVol,quadOrder);
        end

        function ddE = computeDeviatoricEnergyDamageHessian(obj,u,phi,quadOrder)
            ddmu   = obj.d2mu(phi.fun);
            ddeDev = DeviatoricElasticEnergyDensity(u,ddmu);
            ddE    = obj.computeMassWithFunction(ddeDev,quadOrder);
        end

        function Mf = computeMassWithFunction(obj,f,quadOrder)
            Mf = IntegrateLHS(@(u,v) f.*DP(v,u),obj.testPhi,obj.testPhi,obj.mesh,'Domain',quadOrder);
        end

        function V = obtainTensorVolumetric(obj,phi)
            N     = obj.mesh.ndim;
            kappa = Expand(obj.k(phi.fun),4);
            IxI   = ConstantFunction.create(kronEye(N),obj.mesh);
            V     = kappa.*IxI;
        end

        function D = obtainTensorDeviatoric(obj,phi)
            N   = obj.mesh.ndim;
            muE = Expand(obj.mu(phi.fun),4);
            lam = Expand(-(2/N)*obj.mu(phi.fun),4);
            I   = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI = ConstantFunction.create(kronEye(N),obj.mesh);
            D   = 2*muE.*I + lam.*IxI;
        end
        
    end
    
end