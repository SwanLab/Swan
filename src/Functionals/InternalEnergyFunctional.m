classdef InternalEnergyFunctional < handle
    
    properties (Access = private)
        mesh
        material
        testPhi
        testU
    end
    
    methods (Access = public)
        
        function obj = InternalEnergyFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunctional(obj,u,phi,quadOrder)
            C = obj.material.obtainTensor(phi);
            energyFun = DDP(SymGrad(u),DDP(C,SymGrad(u)));
            int  = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);
        end

        function Ju = computeGradientDisplacement(obj,u,phi,quadOrder)  
            C = obj.material.obtainTensor(phi);
            sigma = DDP(C,SymGrad(u));

            s.mesh = obj.mesh;
            s.type = 'ShapeSymmetricDerivative';
            s.quadratureOrder = quadOrder;  
            RHS = RHSIntegrator.create(s);
            Ju  = RHS.compute(sigma,obj.testU);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            dC = obj.material.obtainTensorDerivative(phi);
            dEnergyFun = DDP(SymGrad(u),DDP(dC,SymGrad(u)));
            
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS  = RHSIntegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,obj.testPhi);
        end

        function Huu = computeHessianDisplacement(obj,u,phi,quadOrder)   
            C = obj.material.obtainTensor(phi);
            
            s.trial    = u;
            s.test     = u;
            s.material = C;
            s.mesh     = obj.mesh;
            s.type     = 'ElasticStiffnessMatrix';
            s.quadratureOrder = quadOrder;
            LHS = LHSIntegrator.create(s);
            Huu = LHS.compute();
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)  
            ddC = obj.material.obtainTensorSecondDerivative(phi);
            ddEnergyFun = DDP(SymGrad(u),DDP(ddC,SymGrad(u)));
            
            s.function = ddEnergyFun;
            s.trial    = obj.testPhi;
            s.test     = obj.testPhi;
            s.mesh     = obj.mesh;
            s.type     = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSIntegrator.create(s);
            Hphiphi = 0.5*LHS.compute();
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