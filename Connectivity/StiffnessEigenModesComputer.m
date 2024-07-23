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
            M  = obj.computeMassMatrixWithFunction(m);
            dK = obj.createStiffnessMatrixWithFunction(dalpha);
            dM = obj.computeMassMatrixWithFunction(dm);            
            [lambda,phi] = obj.obtainLowestEigenValuesAndFunction(K,M);
            dlambda = phi'*dK*phi - lambda*phi'*dM*phi;
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
            isRight = @(coor) abs(coor(:,2))==xMax;

            isDir   = @(coor)  isDown(coor) | isUp(coor) | isLeft(coor) | isRight(coor);
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1];
            sDir{1}.value     = 0;

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
            s.function        = obj.createCompositeFunction(fun);
            s.type            = 'StiffnessMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            K = obj.fullToReduced(K);
        end

        function f = createCompositeFunction(obj,fun)
            s.l2function     = obj.density;
            s.handleFunction = fun;
            s.mesh           = obj.mesh;
            f = CompositionFunction(s);
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
            s.function = obj.createCompositeFunction(fun);
            s.quadratureOrder = 2;
            s.type            = 'MassMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();   
            M = obj.fullToReduced(M);
        end       
                
        function [eigV1,eigF1] = obtainLowestEigenValuesAndFunction(obj,K,M)

            [eigF,eigV] = eigs(K,M,4,'smallestabs');
            eigV1 = eigV(1);
            eigF1 = eigF(:,1);
            

        %    [V,eigLHSNewman]  = eigs(K,[],10,'smallestabs');
         %   [Vr,eigLHSDirichlet] = eigs(Kr,[],10,'smallestabs');


            % fV = zeros(size(V(:,1)));
            % fV(obj.boundaryConditions.dirichlet,1) = obj.boundaryConditions.dirichlet_values;
            % fV(obj.boundaryConditions.free,1) = Vr(:,1);
            % s.fValues = fV;
            % s.mesh    = obj.mesh;
            % vV = P1Function(s);
            % vV.plot()
            % 
            % 
            % fV = V(:,2);
            % s.fValues = fV;
            % s.mesh    = obj.mesh;
            % vN = P1Function(s);
            % vN.plot()
    
        end   


        
    end
    
end