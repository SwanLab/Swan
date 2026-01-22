classdef StiffnessEigenModesComputer < handle
    
    properties (Access = public)
        eigsCell
    end
    
    properties (Access = private)
        material 
        massInterpolator
        boundaryConditions
        epsilon
        p
        dim
        phiOld
    end
    
    properties (Access = private)
        mesh
        density
        test
        trial
        phi
    end
    
    methods (Access = public)
        
        function obj = StiffnessEigenModesComputer(cParams)
            obj.init(cParams)                      
        end

        function [lambda,dlambda]  = computeFunctionAndGradient(obj,xR)
            obj.material.setDesignVariable({xR});
            K  = obj.createStiffnessMatrixWithFunction();
            M  = obj.computeMassMatrixWithFunction(xR);
            Kreduced = obj.fullToReduced(K);
            Mreduced = obj.fullToReduced(M);
            [lambda, phi] = obj.obtainLowestEigenValuesAndFunction(Kreduced, Mreduced, 4);          
            dlambda = obj.computeLowestEigenValueGradient(lambda, phi, xR); 
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.test  = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            obj.trial = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            obj.phi =  LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            obj.material = cParams.material;
            obj.massInterpolator = cParams.massInterpolator;
        end  

        function K = createStiffnessMatrixWithFunction(obj)
            C   = obj.material.obtainTensor();
            f = @(u,v) DDP(SymGrad(v),DDP(C,SymGrad(u)));
            K = IntegrateLHS(f,obj.test,obj.trial, obj.mesh,'Domain',2);
        end

        function M = computeMassMatrixWithFunction(obj,xR)          
            weightFun = obj.massInterpolator.fun;
            weight = obj.createDomainFunction(weightFun,xR);               % xR may be chi or 1 - chi
            f = @(u,v) weight.*DP(v,u);
            M = IntegrateLHS(f,obj.test,obj.trial,obj.mesh,'Domain',2);
        end       

        function f = createDomainFunction(obj,fun,xR)
            s.operation = @(xV) obj.createConductivityAsDomainFunction(fun,xR,xV);
            s.mesh      = obj.mesh;
            f = DomainFunction(s);
        end

        function fV = createConductivityAsDomainFunction(obj,fun,xR,xV)
            densV = xR.evaluate(xV);
            fV = fun(densV);
        end

        function K = fullToReduced(obj,K)
            sS.type      = 'DIRECT';
            %sS.type      = 'CG';
            solver       = Solver.create(sS);
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solver     = solver;
            s.boundaryConditions = obj.boundaryConditions;
            s.BCApplier      = obj.createBCApplier();
            ps    = ProblemSolver(s);  
            K = ps.full2Reduced(K);
        end

        function bc = createBCApplier(obj)
            s.mesh = obj.mesh;
            s.boundaryConditions = obj.boundaryConditions;
            bc = BCApplier(s);            
        end
                
        function [eigV1,eigF1] = obtainLowestEigenValuesAndFunction(obj,K,M,n)
            [eigF,eigV] = eigs(K,M,n,'smallestabs');
            eigV1 = eigV(1,1);
            eigF1 = eigF(:,1);
            disp(diag(eigV))
        end   

        function fV = fillVectorWithHomogeneousDirichlet(obj,phi)
            ndofs = obj.phi.nDofs;           
            fV = zeros(ndofs,1);
            dofsDir = obj.boundaryConditions.dirichlet_dofs;
            fV(dofsDir,1) = obj.boundaryConditions.dirichlet_vals;
            free = setdiff(1:ndofs,obj.boundaryConditions.dirichlet_dofs);
            fV(free,1) = phi;
        end

        function dlambda = computeLowestEigenValueGradient(obj, lambda, phi, xR)
            dalpha = obj.createDomainFunction(obj.elasticityInterpolator.dfun, xR);
            dm     = obj.createDomainFunction(obj.massInterpolator.dfun, xR);  
            fValues = obj.fillVectorWithHomogeneousDirichlet(phi);
            obj.phi.setFValues(fValues);
            dlambda = (dalpha.*DP(Grad(obj.phi), Grad(obj.phi)) - lambda*dm.*obj.phi.*obj.phi); 
        end
        
    end
    
end