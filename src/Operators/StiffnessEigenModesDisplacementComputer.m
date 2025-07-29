classdef StiffnessEigenModesDisplacementComputer < handle
    
    properties (Access = public)
        phiDCont
        eigsCell
    end
    
    properties (Access = private)
        conductivity 
        massInterpolator
        boundaryConditions
        epsilon
        p
        phiOld
    end
    
    properties (Access = private)
        mesh
        density
        test
        trial
        eigenF
        material
    end
    
    methods (Access = public)
        
        function obj = StiffnessEigenModesDisplacementComputer(cParams)
            obj.init(cParams)  
            obj.createConductivityInterpolator();   
            obj.createMassInterpolator();                        
        end

        function [eigsV, eigsF] = getEigenModesComputer(obj, dens, n)
            [Kreduced, Mreduced] = createMatrices(obj, dens);
            [eigV, eigF] = obj.obtainEigenValuesAndFunction(Kreduced, Mreduced, n);
            eigsF = [];
            for i = 1:n
                obj.DirichletEigenModeToLagrangianFunction(eigF(:,i));
                phi = obj.eigenF;
                eigsF = [eigsF, phi];
            end
            eigsV = [diag(eigV)];
        end

        function [Kreduced,Mreduced] = createMatrices(obj, dens)
            obj.density = dens;
            alpha  = obj.conductivity.fun;
            m      = obj.massInterpolator.fun;
            K  = obj.createStiffnessMatrixWithFunction(alpha);
            Kreduced = obj.fullToReduced(K);
            M  = obj.computeMassMatrixWithFunction(m);
            Mreduced = obj.fullToReduced(M);
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.test  = LagrangianFunction.create(obj.mesh,3,'P1');
            obj.trial = LagrangianFunction.create(obj.mesh,3,'P1');
            obj.eigenF =  LagrangianFunction.create(obj.mesh,3,'P1'); %posava 1 abans
            obj.material = cParams.material;
        end  
        
        function createConductivityInterpolator(obj)
            s.interpolation  = 'SIMPThermal';   
            s.f0   = 1e-3;                                             
            s.f1   = 1;                                                    
            s.pExp = 2;
            a = MaterialInterpolator.create(s);
            obj.conductivity = a;            
        end            

        function createMassInterpolator(obj)
            s.interpolation  = 'SIMPThermal';                              
            s.f0   = 1e-5;
            s.f1   = 1;
            s.pExp = 1;
            a = MaterialInterpolator.create(s);
            obj.massInterpolator = a;            
        end            

        function K = createStiffnessMatrixWithFunction(obj,fun)
            s.test  = obj.test;
            s.trial = obj.trial;
            s.mesh  = obj.mesh;
            s.material = obj.material;
            s.quadratureOrder = 2;
            s.function        = obj.createDomainFunction(fun);
            s.type            = 'ElasticStiffnessMatrix'; %'StiffnessMatrixWithFunction';
            C   = obj.material.obtainTensor();
            obj.material = C;
            s.material = obj.material;
            lhs = LHSIntegrator.create(s);
            K = lhs.compute();
        end

        function f = createDomainFunction(obj,fun)
            s.operation = @(xV) obj.createConductivityAsDomainFunction(fun,xV);
            s.mesh      = obj.mesh;
            f = DomainFunction(s);
        end

        function fV = createConductivityAsDomainFunction(obj,fun,xV)
            densV = obj.density.evaluate(xV);
            fV = fun(densV);
        end

        function K = fullToReduced(obj,K)
            sS.type      = 'DIRECT';
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

        function M = computeMassMatrixWithFunction(obj,fun)                
            s.test  = obj.test;
            s.trial = obj.trial;
            s.mesh  = obj.mesh;
            s.function = obj.createDomainFunction(fun);
            s.quadratureOrder = 2;
            s.type            = 'MassMatrixWithFunction';
            lhs = LHSIntegrator.create(s);
            M = lhs.compute();   
        end       
        
        function [eigV,eigF] = obtainEigenValuesAndFunction(obj,K,M,n)
            [eigF,eigV] = eigs(K,M,n,'smallestabs');
        end       

        function DirichletEigenModeToLagrangianFunction(obj, eigenF)
            fValues = obj.fillVectorWithHomogeneousDirichlet(eigenF);
            eF      = obj.eigenF;
            fValues = reshape(fValues,[eF.ndimf,eF.nDofs/eF.ndimf])';
            obj.eigenF.setFValues(fValues);
        end

        function fV = fillVectorWithHomogeneousDirichlet(obj,eigenF)
            ndofs = obj.eigenF.nDofs;           
            fV = zeros(ndofs,1);
            dofsDir = obj.boundaryConditions.dirichlet_dofs;
            fV(dofsDir,1) = obj.boundaryConditions.dirichlet_vals;
            free = setdiff(1:ndofs,obj.boundaryConditions.dirichlet_dofs);
            fV(free,1) = eigenF;
        end

    end
    
end