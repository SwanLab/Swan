classdef StiffnessEigenModesComputer < handle
    
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
    end
    
    methods (Access = public)
        
        function obj = StiffnessEigenModesComputer(cParams)
            obj.init(cParams)  
%             obj.createBoundaryConditions();    
            obj.createConductivityInterpolator();   
            obj.createMassInterpolator();                        
        end

        function [lambda,dlambda]  = computeFunctionAndGradient(obj,dens)  
            obj.density = dens;
            alpha  = obj.conductivity.fun;
            dalpha = obj.conductivity.dfun;
            m      = obj.massInterpolator.fun;
            dm     = obj.massInterpolator.dfun;            
            K  = obj.createStiffnessMatrixWithFunction(alpha);
            Kreduced = obj.fullToReduced(K);
            M  = obj.computeMassMatrixWithFunction(m);
            Mreduced = obj.fullToReduced(M);
            [lambdaD, phiD] = obj.obtainLowestEigenValuesAndFunction(Kreduced, Mreduced, 1);
            %[lambdaN, phiD] = obj.obtainLowestEigenValuesAndFunction(K,M,1);
            % obj.plotDirichletEigenMode(phiD)
            % obj.plotNeumannEigenMode(phiN)
            dalphaCont = obj.createDomainFunction(dalpha);
            dmCont     = obj.createDomainFunction(dm);
            obj.DirichletEigenModeToLagrangianFunction(phiD);
            obj.phiDCont = obj.eigenF;
            lambda  = lambdaD;
            dlambda = obj.computeLowestEigenValueGradient(dalphaCont, dmCont, obj.phiDCont, lambdaD); 
        end

        function [eigsV, eigsF] = getEigenModesComputer(obj, dens, n)
            [Kreduced, Mreduced] = createMatrices(obj, dens);
            [eigV, eigF] = obj.obtainEigenValuesAndFunction(Kreduced, Mreduced, n);
            eigsF = [];
            for i = 1:n
%                 fV = eigF(:,i);
%                 s.fValues = fV;
%                 s.mesh    = obj.mesh;
%                 s.order   = 'P1';
%                 phi = LagrangianFunction(s);
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
            obj.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.eigenF =  LagrangianFunction.create(obj.mesh,1,'P1');
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
            s.quadratureOrder = 2;
            s.function        = obj.createDomainFunction(fun);
            s.type            = 'StiffnessMatrixWithFunction';
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
                
        function [eigV1,eigF1] = obtainLowestEigenValuesAndFunction(obj,K,M,n)
            [eigF,eigV] = eigs(K,M,4,'smallestabs');
            i = 1;
            eigV1 = eigV(i,i);
            eigF1 = eigF(:,i);
            if i ~= 1
                disp('SWITCHING, i ='+string(i)+'lambda = '+string(eigV1))
            end
            if abs((eigV(1,1) - eigV(2,2))/eigV(2,2)) < 0.01
                disp('MULTIPLICITY')
            end
            disp(eigV1)
        end   

        function [i] = modalAssuranceCriterion(obj,eigF)
            if ~isempty(obj.phiOld)
                old = obj.phiOld;
                num = (old'*eigF).^2;
                den = dot(old,old)*dot(eigF,eigF);
                mac = num./den;
                [m, i] = max(mac);
                if m < 0.7
                    obj.phiOld = eigF(:,i);
                    disp('tracked eigenmode updated')
                end
            else
                obj.phiOld = eigF(:,1);
                i = 1;
            end
        end   

        function [eigV,eigF] = obtainEigenValuesAndFunction(obj,K,M,n)
            [eigF,eigV] = eigs(K,M,n,'smallestabs');
        end   

        function plotNeumannEigenMode(obj,eigenF)
            fV = eigenF;
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            vV = LagrangianFunction(s);
            vV.plot()
        end        

        function DirichletEigenModeToLagrangianFunction(obj, eigenF)
            fValues = obj.fillVectorWithHomogeneousDirichlet(eigenF);
            obj.eigenF.setFValues(fValues);
        end

        function plotDirichletEigenMode(obj,eigenF)
            obj.DirichletEigenModeToLagrangianFunction(eigenF).plot()
        end

        function fV = fillVectorWithHomogeneousDirichlet(obj,eigenF)
            ndofs = obj.eigenF.nDofs;           
            fV = zeros(ndofs,1);
            dofsDir = obj.boundaryConditions.dirichlet_dofs;
            fV(dofsDir,1) = obj.boundaryConditions.dirichlet_vals;
            free = setdiff(1:ndofs,obj.boundaryConditions.dirichlet_dofs);
            fV(free,1) = eigenF;
        end

        function dlambda = computeLowestEigenValueGradient(obj, dalpha, dm, phi, lambda)
            dlambda = (dalpha.*DP(Grad(phi), Grad(phi)) - lambda*dm.*phi.*phi); 
        end
        

    end
    
end