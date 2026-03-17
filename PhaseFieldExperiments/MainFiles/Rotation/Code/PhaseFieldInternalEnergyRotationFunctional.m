classdef PhaseFieldInternalEnergyRotationFunctional < handle
    
    properties (Access = private)
        mesh
        material
        testPhi
        testU
    end
    
    methods (Access = public)
        
        function obj = PhaseFieldInternalEnergyRotationFunctional(cParams)
            obj.init(cParams)
        end

        function F = computeCost(obj,u,theta,phi,quadOrder)
            C = obj.material.obtainTensor(phi);
            eps = SymGrad(u); epsRot = obj.rotateTensor(eps,theta);
            sig = DDP(C,epsRot); sigRot = obj.rotateTensor(sig,-theta);
            energyFun = DDP(eps,sigRot);
            int  = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);
        end

        function Ju = computeGradientDisplacement(obj,u,theta,phi,quadOrder)  
            C = obj.material.obtainTensor(phi);
            eps = SymGrad(u); epsRot = obj.rotateTensor(eps,theta);
            sig = DDP(C,epsRot); sigRot = obj.rotateTensor(sig,-theta);
            Ju = IntegrateRHS(@(v) DDP(SymGrad(v),sigRot),obj.testU,obj.mesh,'Domain',quadOrder);
        end

        function Jphi = computeGradientDamage(obj,u,theta,phi,quadOrder)
            dC = obj.material.obtainTensorDerivative(phi);
            eps = SymGrad(u); epsRot = obj.rotateTensor(eps,theta);
            sig = DDP(dC,epsRot); sigRot = obj.rotateTensor(sig,-theta);
            dEnergyFun = DDP(eps,sigRot);
            Jphi = IntegrateRHS(@(v) (1/2)*DP(v,dEnergyFun),obj.testPhi,obj.mesh,'Domain',quadOrder);
        end

        function Huu = computeHessianDisplacement(obj,~,theta,phi,quadOrder)   
            C = obj.material.obtainTensor(phi);
            eps = @(u) SymGrad(u); epsRot = @(u) obj.rotateTensor(eps(u),theta);
            sig = @(u) DDP(C,epsRot(u)); sigRot = @(u) obj.rotateTensor(sig(u),-theta);
            Huu = IntegrateLHS(@(u,v) DDP(SymGrad(v),sigRot(u)),obj.testU,obj.testU,obj.mesh,'Domain',quadOrder);
        end

        function Hphiphi = computeHessianDamage(obj,u,theta,phi,quadOrder)  
            ddC = obj.material.obtainTensorSecondDerivative(phi);
            eps = SymGrad(u); epsRot = obj.rotateTensor(eps,theta);
            sig = DDP(ddC,epsRot); sigRot = obj.rotateTensor(sig,-theta);
            ddEnergyFun = DDP(eps,sigRot);
            Hphiphi = IntegrateLHS(@(u,v) (1/2)*ddEnergyFun.*DP(u,v),obj.testPhi,obj.testPhi,obj.mesh,'Domain',quadOrder);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;            
            obj.testPhi  = copy(cParams.testSpace.phi);
            obj.testU    = copy(cParams.testSpace.u);
        end
        
        function tensorRot = rotateTensor(obj,tensor,angle)
            Rot = obj.computeRotationMatrix(angle);
            tensorRot = Rot*tensor*Rot';
        end

        function Rot = computeRotationMatrix(obj,angle)
            s.operation = @(xV) [cos(Expand(angle,2).evaluate(xV)), -sin(Expand(angle,2).evaluate(xV));
                                 sin(Expand(angle,2).evaluate(xV))   cos(Expand(angle,2).evaluate(xV))];
            s.ndimf = 4;
            s.mesh = obj.mesh;
            Rot = DomainFunction(s);
        end
    end
    
end