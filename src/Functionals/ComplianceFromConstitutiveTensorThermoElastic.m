classdef ComplianceFromConstitutiveTensorThermoElastic < handle

    properties (Access = private)
        quadrature
    end

    properties (Access = private)
        mesh
        stateProblem
        adjointProblem
        materialInterpolator
        T0      
    end

    methods (Access = public)
        function obj = ComplianceFromConstitutiveTensorThermoElastic(cParams)
            obj.init(cParams);
            obj.createQuadrature();
            obj.createAdjointProblem();
            obj.createAdjointBoundaryConditions();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,C,dC,kappa, dkappa)
            [u,T] = obj.computeStateVariable(C,kappa);
            [p] = obj.computeAdjointVariable(kappa,C,u);
            J  = obj.computeFunction(C,u);
            dJ = obj.computeGradient(dC,u,dkappa,T,p);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.stateProblem = cParams.stateProblem;
            obj.T0           = cParams.T0;
        end

        function createQuadrature(obj)
            quad = Quadrature.create(obj.mesh,2);
            obj.quadrature = quad;
        end

        function createAdjointProblem(obj)
            s.mesh = obj.mesh;
            s.conductivity = obj.materialInterpolator; 
            Q = LagrangianFunction.create(obj.mesh,1,'P1');
            fValues = ones(Q.nDofs,1);
            Q.setFValues(fValues);
            s.source       = Q;  
            s.dim = '2D';
            s.boundaryConditions = obj.createAdjointBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';         
            obj.adjointProblem = ThermalProblem(s); 
            
        end

        function [u,T] = computeStateVariable(obj,C,kappa)
            obj.stateProblem.updateMaterial(C,kappa);
            obj.stateProblem.solve();
            u = obj.stateProblem.uFun;
            T = obj.stateProblem.tFun;
        end

        function [p] = computeAdjointVariable(obj,kappa,C,u)
            I = ConstantFunction.create(eye(2),obj.mesh);
            beta = obj.stateProblem.alpha.*DDP(C,I);    
            newSource = -DDP(beta,SymGrad(u));        % updating RHS
            obj.adjointProblem.updateSource(newSource);
            obj.adjointProblem.solve(kappa);
            p = obj.adjointProblem.uFun;
        end

        function bcAdj = createAdjointBoundaryConditions(obj)
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            isDir   = @(coor) abs(coor(:,2))==yMin & abs(coor(:,1))>=0.4*xMax & abs(coor(:,1))<=0.6*xMax;  
            
            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = 1;
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
            bcAdj = BoundaryConditions(s); 
        end

        function J = computeFunction(obj,C,u)
            dCompliance = ElasticEnergyDensity(C,u);
            J           = Integrator.compute(dCompliance,obj.mesh,obj.quadrature.order);
        end

        function dj = computeGradient(obj,dC,u,dkappa,T,p)
            nDesVar = length(dC);
            dj      = cell(nDesVar,1);
            for i = 1:nDesVar
                strain  = SymGrad(u);
                dStress = DDP(dC{i},strain);
                I = ConstantFunction.create(eye(2),obj.mesh);
                dbeta = obj.stateProblem.alpha.*DDP(dC{i},I);
                dj{i}   = -0.5.*DDP(strain, dStress) + DDP(dbeta.*(T-obj.T0),SymGrad(u)) + DP(Grad(T),dkappa.*Grad(p)); 
            end
        end

    end
end