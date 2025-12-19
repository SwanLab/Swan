classdef GripperLevelSetCircle < handle

    properties (Access = private)
        filename
        mesh
        filter
        filterConnectivity
        filterPerimeter
        filterAdjointConnectivity 
        designVariable
        materialInterpolator
        conductivityInterpolator
        massInterpolator
        physicalProblem
        compliance
        perimeter
        volume
        cost
        constraint
        primalUpdater
        optimizer
        minimumEigenValue
        eigenvalue
    end

    methods (Access = public)

        function obj = GripperLevelSetCircle()
            obj.init()
            obj.createMesh();
            obj.createDesignVariable();
            obj.createFilter();
            obj.createFilterConnectivity();
            obj.createFilterPerimeter();
            obj.createMaterialInterpolator();
            obj.createElasticProblem();
            obj.createCompliance();
            obj.createConductivityInterpolator();
            obj.createMassInterpolator();
            obj.createPerimeter();
            obj.createEigenValueConstraint();
            obj.createEigenValue();
            obj.createVolumeConstraint();
            obj.createCost();
            obj.createConstraint();
            obj.createPrimalUpdater();
            obj.createOptimizer();
            saveas(figure(1),'GripperLevelSetCircle.fig');
            saveas(figure(2),'GripperLevelSetCircleMonitoring.fig');
            obj.designVariable.fun.print('GripperLevelSetCirclefValues');
            writematrix(obj.designVariable.fun.fValues,'GripperLevelSetCirclefValues.txt');
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'Gripping';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function createDesignVariable(obj)
%             s.type = 'Full';
%             g      = GeometricalFunction(s);
            lsFun = LagrangianFunction.create(obj.mesh,1,'P1');
            lsFun.setFValues(importdata('GripperLevelSetCirclefValues.txt'))
%             lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            h = obj.mesh.computeMinCellSize();
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
%             f.updateEpsilon(1.0*h);
            obj.filter = f;
        end

       function createFilterPerimeter(obj)
            s.filterType = 'PDE';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            f.updateEpsilon(1*obj.mesh.computeMinCellSize());
            obj.filterPerimeter = f;           
        end



        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filterConnectivity = f;
% 
% %             s.filterType = 'FilterAndProject';
%             s.filterType = 'CloseOperator';   
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'PDE';
%             s.beta       = 10.0; % 4 2 
%             s.eta        = 0.2;
%             obj.filterConnectivity = Filter.create(s);
            
%             s.filterType = 'FilterAdjointAndProject';

%             s.filterType = 'CloseAdjointOperator';   
%             s.mesh       = obj.mesh;
%             s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
%             s.filterStep = 'PDE';
%             s.beta       = 10.0; 
%             s.eta        = 0.2;
%             obj.filterAdjointConnectivity = Filter.create(s);
        end

        function createMaterialInterpolator(obj)
            E0 = 1e-3;
            nu0 = 1/3;
            ndim = obj.mesh.ndim;
            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);


            E1 = 1;
            nu1 = 1/3;
            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)
            f = obj.designVariable.fun;           
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '2D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = DirectSolver();
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function createCompliance(obj)
            s.mesh         = obj.mesh;
            s.filter       = obj.filter;
            s.material     = obj.createMaterial();
            s.stateProblem = obj.physicalProblem;
            s.filename     = obj.filename;
            c = NonSelfAdjointComplianceFunctional(s);
            obj.compliance = c;
        end

        function createConductivityInterpolator(obj) 
            s.interpolation  = 'SimpAllThermal';
            s.f0   = 1e-3;                                             
            s.f1   = 1;  
            s.dim  = '2D';
            a = MaterialInterpolator.create(s);
            obj.conductivityInterpolator = a;            
        end 

        function createMassInterpolator(obj)
            s.interpolation  = 'SIMPThermal';                              
            s.f0   = 1e-3;
            s.f1   = 1;
            s.pExp = 1;
            a = MaterialInterpolator.create(s);
            obj.massInterpolator = a;            
        end   

        function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnectivity;
            s.filterAdjoint     = obj.filterAdjointConnectivity;
            s.BC = 'Neumann';
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.conductivityInterpolator = obj.conductivityInterpolator; 
            s.massInterpolator         = obj.massInterpolator; 
            s.targetEigenValue  = 6.0;     
            s.isCompl           = false;
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
        end

        function createEigenValue(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnectivity;     
            s.filterAdjoint     = obj.filterAdjointConnectivity;
            s.BC = 'Neumann';
            s.conductivityInterpolator = obj.conductivityInterpolator; 
            s.massInterpolator         = obj.massInterpolator; 
            s.isCompl           = false;
            s.dim               = '2D';
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            obj.eigenvalue = MaximumEigenValueFunctional(s);
        end

        function uMesh = createBaseDomain(obj)
            levelSet         = -ones(obj.mesh.nnodes,1);
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            uMesh = UnfittedMesh(s);
            uMesh.compute(levelSet);
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.test = LagrangianFunction.create(obj.mesh,1,'P1');
            s.volumeTarget = 0.6; 
            s.uMesh = obj.createBaseDomain();
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createPerimeter(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterPerimeter;
            p = SimplePerimeterFunctional(s);
            obj.perimeter = p;
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.shapeFunctions{2} = obj.eigenvalue;
            s.weights           = [1.0,1.0]; 
            s.Msmooth           = obj.createMassMatrix();
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
%             s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.maxIter        = 4000;
            s.tolerance      = 1e-8;
%             s.constraintCase = {'INEQUALITY','INEQUALITY'};
            s.constraintCase = {'INEQUALITY'};
            s.etaNorm        = 0.05;
            s.etaNormMin = 0.05;
            s.etaMax = 10.0;
            s.etaMaxMin = 0.1;
            s.gJFlowRatio    = 0.01;
%             s.gJFlowRatio = 0.1;
            s.primalUpdater  = obj.primalUpdater;
            s.gif = true;
            s.gifName = char("gripper");
            s.printing = false;
            s.printName = [];
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
        end

        function bc = createBoundaryConditions(obj)
              femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            sPL       = obj.computeCondition(s.pointload);
            sDir      = obj.computeCondition(s.dirichlet);

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
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function  bc = createEigenvalueBoundaryConditions(obj)
            xMin    = min(obj.mesh.coord(:,1));
            yMin    = min(obj.mesh.coord(:,2));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isFront = @(coor) abs(coor(:,2))==yMin;
            isBack = @(coor) abs(coor(:,2))== yMax;

            isTopGripper    = @(coor)  (abs(coor(:,1)) < 0.1  & coor(:,2) == 0.6 );
            isBottomGripper = @(coor)  (abs(coor(:,1)) < 0.1  & coor(:,2) == 0.4 );
            isRightGripper = @(coor)  (abs(coor(:,1)) - 0.1 < 1e-6 & coor(:,2) > 0.4 & coor(:,2) < 0.6);

            isTopRight      = @(coor)  (abs(coor(:,1)) >= 0.92 & coor(:,2) == 1 );
            isBottomRight   = @(coor)  (abs(coor(:,1)) >= 0.92 & abs(coor(:,2)) <= 1e-8 ); % not exactly 0 in the mesh

            isDir   = @(coor) isLeft(coor) | isRight(coor) | isFront(coor) | isBack(coor) | isTopGripper(coor) | isBottomGripper(coor) | isRightGripper(coor);  
%             isDir   = @(coor) isTopGripper(coor) | isBottomGripper(coor) | isRightGripper(coor);  
% %             isDir   = @(coor) isTopRight(coor) | isBottomRight(coor);  
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
            bc = BoundaryConditions(s);  
        end
    end


    methods (Static, Access=private)
        function sCond = computeCondition(conditions)
            nodes = @(coor) 1:size(coor,1);
            dirs  = unique(conditions(:,2));
            j     = 0;
            for k = 1:length(dirs)
                rowsDirk = ismember(conditions(:,2),dirs(k));
                u        = unique(conditions(rowsDirk,3));
                for i = 1:length(u)
                    rows   = conditions(:,3)==u(i) & rowsDirk;
                    isCond = @(coor) ismember(nodes(coor),conditions(rows,1));
                    j      = j+1;
                    sCond{j}.domain    = @(coor) isCond(coor);
                    sCond{j}.direction = dirs(k);
                    sCond{j}.value     = u(i);
                end
            end
        end

    end

end