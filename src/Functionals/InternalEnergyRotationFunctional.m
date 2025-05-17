classdef InternalEnergyRotationFunctional < handle
    
    properties (Access = private)
        mesh
        material
        testPhi
        testU
    end
    
    methods (Access = public)
        
        function obj = InternalEnergyRotationFunctional(cParams)
            obj.init(cParams)            
        end
        
        function F = computeFunction(obj,u,phi,quadOrder)
            alpha = obj.updateRotation(u,phi);
            C = obj.material.obtainTensor(phi);
            sigma = DDP(C{1},SymGrad(u));
            energyFun = DDP(SymGrad(u),sigma);
            int = Integrator.create('Function',obj.mesh,quadOrder);
            F = 0.5*int.compute(energyFun);

            obj.computeStressPrincipalDirections(sigma)
        end

        function computeStressPrincipalDirections(obj,sigma)

        end

        function Ju = computeGradientDisplacement(obj,u,phi,quadOrder)
            alpha = obj.updateRotation(u,phi);
            C = obj.material.obtainTensor(phi);
            sigma = DDP(C{1},SymGrad(u));
            test = obj.testU;

            s.mesh = obj.mesh;
            s.quadratureOrder = quadOrder;
            s.type = 'ShapeSymmetricDerivative';
            RHS = RHSIntegrator.create(s);
            Ju = RHS.compute(sigma,test);
        end

        function Jphi = computeGradientDamage(obj,u,phi,quadOrder)
            alpha = obj.updateRotation(u,phi);
            dC = obj.material.obtainTensorDerivative(phi);
            dEnergyFun = DDP(SymGrad(u),DDP(dC{1},SymGrad(u)));
            test = obj.testPhi;
            
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            RHS = RHSIntegrator.create(s);
            Jphi = 0.5*RHS.compute(dEnergyFun,test);
        end

        function Huu = computeHessianDisplacement(obj,u,phi,quadOrder)
            alpha = obj.updateRotation(u,phi);
            C = obj.material.obtainTensor(phi); 
            s.type     = 'ElasticStiffnessMatrix';
            s.mesh     = obj.mesh;
            s.trial     = u;
            s.test      = u;
            s.material = C{1};
            s.quadratureOrder = quadOrder;
            LHS = LHSIntegrator.create(s);
            Huu = LHS.compute();
        end

        function Hphiphi = computeHessianDamage(obj,u,phi,quadOrder)
            alpha = obj.updateRotation(u,phi);
            ddC = obj.material.obtainTensorSecondDerivative(phi);
            ddEnergyFun = DDP(SymGrad(u),DDP(ddC{1},SymGrad(u)));
            
            s.fun = ddEnergyFun;
            s.trial = obj.testPhi;
            s.test  = obj.testPhi;
            s.mesh = obj.mesh;
            s.type = 'MassMatrixWithFunction';
            s.quadratureOrder = quadOrder;
            LHS = LHSIntegrator.create(s);
            Hphiphi = 0.5*LHS.compute();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.material = cParams.material;            
            obj.testPhi = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.testU   = LagrangianFunction.create(obj.mesh, 2, 'P1');
        end
        
        function alpha = updateRotation(obj,u,phi)
            s.type = '2D';
            s.ndim = 2;
            s.eigenValueComputer.type = 'PRECOMPUTED';
            p = PrincipalDirectionComputer.create(s);
            
            
            %sigTensor = permute(sigTensor,[2 1 3]);
            p.compute([0 0 1])
            %p.compute(rand(3,4,10))
            alpha = p.direction;
        end
    end
    
end