classdef MultimaterialTesting < handle

    properties (Access = private)
        E
        nu
        mesh
        designVariable
        filter
        boundaryConditions
        simpAlls
        physicalProblem
        compliance
        volumeA
        volumeB
        volumeC
        cost
        constraint
        dualVariable
        primalUpdater
        optimizer
    end

    methods (Access = public)

        function obj = MultimaterialTesting()
            obj.init()
            obj.createMesh();
            obj.createFilter();
            obj.createDesignVariable();
            obj.createBoundaryConditions();
            obj.createElasticProblem();
            obj.createSimpAlls();
            obj.createCompliance();
            obj.createVolumeConstraints();
            obj.createCost();
            obj.createConstraint();
            obj.createDualVariable();
            obj.createPrimalUpdater();
            obj.createOptimizer();
        end

        function lsVals = solve(obj)
            obj.optimizer.solveProblem();
            lsVals = obj.designVariable.fun.fValues;
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
            obj.E    = [1,0.5,0.25,1e-3];
            obj.nu   = (1/3)*[1,1,1,1];
        end

        function createMesh(obj)
            obj.mesh = TriangleMesh(2,1,14,7);
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filter   = f;
        end

        function createDesignVariable(obj)
            lsFun{1} = obj.createLevelSetFunction(@(x) -ones(size(x(1,:,:))));
            lsFun{2} = obj.createLevelSetFunction(@(x) -cos(x(1,:,:))+0.5);
            lsFun{3} = obj.createLevelSetFunction(@(x) sin(x(1,:,:))-0.5);
            
            s.type             = 'MultiLevelSet';
            s.lsFun            = lsFun;
            s.mesh             = obj.mesh;
            s.plotting         = true;
            obj.designVariable = DesignVariable.create(s);
        end

        function lsFun = createLevelSetFunction(obj,fH)
            s.type    = 'Given';
            s.fHandle = fH;
            g         = GeometricalFunction(s);
            lsFun     = g.computeLevelSetFunction(obj.mesh);
        end

        function createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.35*yMax & abs(coor(:,2))<=0.65*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'DIRAC');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            obj.boundaryConditions = BoundaryConditions(s);
        end

        function createSimpAlls(obj)
            [muRef,kRef] = obj.computeReferenceShearBulk();
            N = obj.mesh.ndim;
            for i = 1:length(obj.E)-1
                for j = 2:length(obj.E)
                    obj.simpAlls{i,j}.mu  = @(rho) SimpAllInterpolator.computeMu(muRef{i},muRef{j},kRef{i},kRef{j},rho,N);
                    obj.simpAlls{i,j}.k   = @(rho) SimpAllInterpolator.computeKappa(muRef{i},muRef{j},kRef{i},kRef{j},rho,N);
                    obj.simpAlls{i,j}.dmu = @(rho) SimpAllInterpolator.computeMuDerivative(muRef{i},muRef{j},kRef{i},kRef{j},rho,N);
                    obj.simpAlls{i,j}.dk  = @(rho) SimpAllInterpolator.computeKappaDerivative(muRef{i},muRef{j},kRef{i},kRef{j},rho,N);
                end
            end
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = [];
            s.dim = '2D';
            s.boundaryConditions = obj.boundaryConditions;
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstiutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filter;
            s.complainceFromConstitutive = obj.createComplianceFromConstiutive();
            s.C                          = obj.computeMaterial();
            s.dC                         = obj.computeMaterialDerivative();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function C = computeMaterial(obj)
            [muRef,kRef] = obj.computeReferenceShearBulk();
            N            = obj.mesh.ndim;
            [mu,kappa]   = MultiMaterialInterpolator.computeMuKappa(muRef,kRef,obj.simpAlls,N);
            lambda       = @(x) LameParametersConverter.computeLambdaFromBulkAndShear(kappa(x),mu(x),N);
            mu           = @(x) Expand(mu(x),4);
            lambda       = @(x) Expand(lambda(x),4);
            I            = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI          = ConstantFunction.create(kronEye(N),obj.mesh);
            C            = @(x) 2*mu(x).*I + lambda(x).*IxI;
        end

        function [mu,k] = computeReferenceShearBulk(obj)
            N     = obj.mesh.ndim;
            Evec  = obj.E;
            nuvec = obj.nu;
            mu    = cell(length(Evec),1);
            k     = cell(length(Evec),1);
            for i = 1:length(Evec)
                mu{i} = LameParametersConverter.computeShearFromYoungAndPoisson(Evec(i),nuvec(i));
                k{i}  = LameParametersConverter.computeBulkFromYoungAndPoisson(Evec(i),nuvec(i),N);
            end
        end

        function dC = computeMaterialDerivative(obj)
            N    = obj.mesh.ndim;
            Z    = ConstantFunction.create(0,obj.mesh);
            I    = ConstantFunction.create(1,obj.mesh);
            II   = ConstantFunction.create(eye4D(N),obj.mesh);
            IxI  = ConstantFunction.create(kronEye(N),obj.mesh);
            nMat = length(obj.E);
            dCLoc   = cell(nMat,nMat);
            for i = 1:nMat
                for j = 1:nMat
                    if i==j
                        dmu    = Z;
                        dkappa = Z;
                    elseif i<j
                        dmu    = obj.simpAlls{i,j}.dmu(Z);
                        dkappa = obj.simpAlls{i,j}.dk(Z);
                    else
                        dmu    = - obj.simpAlls{j,i}.dmu(I);
                        dkappa = - obj.simpAlls{j,i}.dk(I);
                    end
                    dlambda = dkappa - (2/N)*dmu;
                    dmu = Expand(dmu,4); dlambda = Expand(dlambda,4);
                    dCLoc{i,j} = 2*dmu.*II + dlambda.*IxI;
                end
            end
            dC = @(x) dCLoc;
        end

        function createVolumeConstraints(obj)
            obj.volumeA = obj.createIndivVolumeConstraint(0.1,1);
            obj.volumeB = obj.createIndivVolumeConstraint(0.1,2);
            obj.volumeC = obj.createIndivVolumeConstraint(0.1,3);
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function v = createIndivVolumeConstraint(obj,target,ID)
            s.volumeTarget = target;
            s.nMat         = 4;
            s.matID        = ID;
            s.mesh         = obj.mesh;
            s.test         = LagrangianFunction.create(obj.mesh,1,'P1');
            s.uMesh        = obj.createBaseDomain();
            v              = MultiMaterialVolumeConstraint(s);
         end

         function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = 1;
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
         end

         function createConstraint(obj)
            s.shapeFunctions{1} = obj.volumeA;
            s.shapeFunctions{2} = obj.volumeB;
            s.shapeFunctions{3} = obj.volumeC;
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
         end

         function createDualVariable(obj)
            s.nConstraints   = 3;
            l                = DualVariable(s);
            obj.dualVariable = l;
         end

         function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

         function createOptimizer(obj)
            s.monitoring     = false;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 3;
            s.tolerance      = 1e-8;
            s.constraintCase = repmat({'EQUALITY'},[3,1]);
            s.primalUpdater  = obj.primalUpdater;
            s.delta          = inf;
            s.etaStar        = 0.02;
            s.gif            = false;
            s.gifName        = [];
            s.printing       = false;
            s.printName      = [];
            obj.optimizer    = OptimizerNullSpace(s);
         end

        function M = createMassMatrix(obj)
            nnodes  = obj.mesh.nnodes*3;
            indices = transpose(1:nnodes);
            vals    = ones(size(indices));
            h       = obj.mesh.computeMeanCellSize();
            M       = h^2*sparse(indices,indices,vals,nnodes,nnodes);
        end


    end
end