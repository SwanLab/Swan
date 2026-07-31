classdef PhaseFieldInternalEnergyFunctional < handle
    
    properties (Access = private)
        mesh
        C
        dC
        d2C
        testPhi
        testU
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldInternalEnergyFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeCost(obj,u,phi,quadOrder)
            Cph = obj.C(phi.fun);
            energyFun = DDP(SymGrad(u),DDP(Cph,SymGrad(u)));
            int  = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);
        end

        function Ju = computeGradientDisplacement(obj,u,phi,quadOrder)  
            Cph = obj.C(phi.fun);
            sigma = DDP(Cph,SymGrad(u));
            Ju = IntegrateRHS(@(v) DDP(SymGrad(v),sigma),obj.testU,obj.mesh,'Domain',quadOrder);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            dCph = obj.dC(phi.fun);
            dEnergyFun = DDP(SymGrad(u),DDP(dCph,SymGrad(u)));
            Jphi = IntegrateRHS(@(v) (1/2)*DP(v,dEnergyFun),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end

        function Huu = computeHessianDisplacement(obj,~,phi,quadOrder)   
            Cph = obj.C(phi.fun);
            Huu = IntegrateLHS(@(u,v) DDP(SymGrad(v),DDP(Cph,SymGrad(u))),obj.testU,obj.testU,obj.mesh,'Domain',quadOrder);
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)  
            d2Cph = obj.d2C(phi.fun);
            ddEnergyFun = DDP(SymGrad(u),DDP(d2Cph,SymGrad(u)));
            Hphiphi = IntegrateLHS(@(u,v) (1/2)*ddEnergyFun.*DP(u,v),obj.testPhi,obj.testPhi,obj.mesh,'Domain',quadOrder);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.C       = cParams.C;
            obj.dC      = cParams.dC;    
            obj.d2C     = cParams.d2C;    
            obj.testPhi = copy(cParams.testSpace.phi);
            obj.testU   = copy(cParams.testSpace.u);
        end
        
    end
    
end