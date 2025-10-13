classdef PhaseFieldInternalEnergyFunctional < handle
    
    properties (Access = private)
        mesh
        material
        testPhi
        testU
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldInternalEnergyFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeCost(obj,u,phi,quadOrder)
            C = obj.material.obtainTensor(phi);
            energyFun = DDP(SymGrad(u),DDP(C,SymGrad(u)));
            int  = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);
        end

        function Ju = computeGradientDisplacement(obj,u,phi,quadOrder)  
            C = obj.material.obtainTensor(phi);
            sigma = DDP(C,SymGrad(u));
            Ju = IntegrateRHS(@(v) DDP(SymGrad(v),sigma),obj.testU,obj.mesh,quadOrder);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            dC = obj.material.obtainTensorDerivative(phi);
            dEnergyFun = DDP(SymGrad(u),DDP(dC,SymGrad(u)));
            Jphi = IntegrateRHS(@(v) (1/2)*DP(v,dEnergyFun),obj.testPhi,obj.mesh,quadOrder);
        end

        function Huu = computeHessianDisplacement(obj,~,phi,quadOrder)   
            C = obj.material.obtainTensor(phi);
            Huu = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u))),obj.testU,obj.testU,obj.mesh,quadOrder);
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)  
            ddC = obj.material.obtainTensorSecondDerivative(phi);
            ddEnergyFun = DDP(SymGrad(u),DDP(ddC,SymGrad(u)));
            Hphiphi = IntegrateLHS(@(u,v) (1/2)*ddEnergyFun.*DP(u,v),obj.testPhi,obj.testPhi,obj.mesh,quadOrder);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;            
            obj.testPhi  = copy(cParams.testSpace.phi);
            obj.testU    = copy(cParams.testSpace.u);
        end
        
    end
    
end