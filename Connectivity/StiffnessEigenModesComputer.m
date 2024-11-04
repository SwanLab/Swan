classdef StiffnessEigenModesComputer < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        conductivity 
        massInterpolator
        boundaryConditions
        Kmatrix
        Mmatrix        
    end
    
    properties (Access = private)
        mesh
        density
    end
    
    methods (Access = public)
        
        function obj = StiffnessEigenModesComputer(cParams)
            obj.init(cParams)  
            obj.createBoundaryConditions();    
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
            dK = obj.createStiffnessMatrixWithFunction(dalpha);
            dM = obj.computeMassMatrixWithFunction(dm);
            [lambdaD,phiD] = obj.obtainLowestEigenValuesAndFunction(Kreduced, Mreduced, 1);
            [lambdaN, phiN] = obj.obtainLowestEigenValuesAndFunction(K,M,1);
            obj.plotDirichletEigenMode(phiD)
            obj.plotNeumannEigenMode(phiN)
            dlambda = dK*phiN - lambdaD*dM*phiN; %phi'*dK*phi - lambda*phi'*dM*phi;
            lambda  = lambdaD;
 
            obj.plotEigenValueDerivative(dlambda)

            % Derivative Expression with PhiD
            phiD_filled = obj.fillVectorWithHomogeneousDirichlet(phiD);
            dlambdaD = dK*phiD_filled - lambdaD*dM*phiD_filled;

            obj.plotEigenValueDerivative(dlambdaD)

            % Derivative Expression from Continuous Version
            dalphaCont = obj.createDomainFunction(dalpha)
            s.fValues = obj.fillVectorWithHomogeneousDirichlet(phiD);
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            phiDCont = LagrangianFunction(s);
            dlambdaCont = obj.computeLowestEigenValueGradient(dalphaCont, phiDCont, lambdaD); 
            dlambdaCont.project('P1',obj.mesh).plot()
        end
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh    = cParams.mesh;
        end

        function createBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDown  = @(coor) abs(coor(:,2))==yMin;
            isUp    = @(coor) abs(coor(:,2))==yMax;
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;

            isDir   = @(coor)  isDown(coor) | isUp(coor) | isLeft(coor) | isRight(coor);  
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1];
            sDir{1}.value     = 0;
            sDir{1}.ndim = 1;

             dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);  
            obj.boundaryConditions = bc;
        end        
        
        function createConductivityInterpolator(obj)
            s.interpolation  = 'SIMPThermal';                              
            s.f0   = 1e-5;                                                 
            s.f1   = 1;                                                    
            s.pExp = 8;
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
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.quadratureOrder = 2;
            s.function        = obj.createDomainFunction(fun);
            s.type            = 'StiffnessMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            % K = obj.fullToReduced(K);
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
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            s.mesh  = obj.mesh;
            s.function = obj.createDomainFunction(fun);
            s.quadratureOrder = 2;
            s.type            = 'MassMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();   
            % M = obj.fullToReduced(M);
        end       
                
        function [eigV1,eigF1] = obtainLowestEigenValuesAndFunction(obj,K,M,n)
            [eigF,eigV] = eigs(K,M,4,'smallestabs');
            eigV1 = eigV(n,n);
            eigF1 = eigF(:,n);
        end   

        function plotNeumannEigenMode(obj,eigenF)
            fV = eigenF;
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            vV = LagrangianFunction(s);
            vV.plot()
        end        

        function plotDirichletEigenMode(obj,eigenF)
%             t = LagrangianFunction.create(obj.mesh,1,'P1');
%             ndofs = t.nDofs;           
%             fV = zeros(ndofs,1);
%             dofsDir = obj.boundaryConditions.dirichlet_dofs;
%             fV(dofsDir,1) = obj.boundaryConditions.dirichlet_vals;
%             free = setdiff(1:ndofs,obj.boundaryConditions.dirichlet_dofs);
%             fV(free,1) = eigenF;
            s.fValues = obj.fillVectorWithHomogeneousDirichlet(eigenF);
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            vV = LagrangianFunction(s);
            vV.plot()
        end

        function fV = fillVectorWithHomogeneousDirichlet(obj,eigenF)
            t = LagrangianFunction.create(obj.mesh,1,'P1');
            ndofs = t.nDofs;           
            fV = zeros(ndofs,1);
            dofsDir = obj.boundaryConditions.dirichlet_dofs;
            fV(dofsDir,1) = obj.boundaryConditions.dirichlet_vals;
            free = setdiff(1:ndofs,obj.boundaryConditions.dirichlet_dofs);
            fV(free,1) = eigenF;
        end

        function plotEigenValueDerivative(obj,dlambda)
            fV = dlambda;
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            vV = LagrangianFunction(s);
            vV.plot()
        end

        function dlambda = computeLowestEigenValueGradient(obj, dalpha, phi, lambda)
            dlambda = dalpha.*DDP(Grad(phi), Grad(phi)) - lambda*phi*phi; % obj.density.* on the second term?
        end
        


    end
    
end