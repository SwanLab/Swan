classdef TopOptLevelSet3DConnectivity< handle

    properties (Access = private)
        mesh
        filter
        filterComp
        filterAdjointComp
        filterConnect
        filterAdjointConnect
        designVariable
        materialInterpolator
        physicalProblem
        compliance
        volume
        cost
        constraint
        dualVariable
        optimizer
        minimumEigenValue
        gJ
        lambda1min
        filterAdjointProj 
        filterProj
        complianceProj
        volumeProj
        eta
        dofsNonDesign
        type
        primalUpdater
    end  

    methods (Access = public)
        function obj = TopOptLevelSet3DConnectivity()
            for type = ["torqueBeam"] %,] %, "cornersSupport","cornersSupport","centralSupport",
                for lambda1min = [2.0]
                    obj.lambda1min = lambda1min;
                    obj.type = convertStringsToChars(type);
                    obj.init()
                    obj.createMesh();
                    obj.createDesignVariable();
                    obj.createFilter();
                    obj.createFilterCompliance();
                    obj.createFilterConnectivity();
                    obj.createMaterialInterpolator();
                    obj.createElasticProblem();
                    obj.createComplianceFromConstitutive();
                    obj.createNonDesignableDomain();
                    obj.createCompliance();
                    obj.createEigenValueConstraint();                             
                    obj.createVolumeConstraint();
                    obj.createCost();
                    obj.createConstraint();
                    obj.createPrimalUpdater();
                    obj.createOptimizer();
                end
            end
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            if isequal(obj.type, 'cornersSupport') 
                obj.mesh = HexaMesh(1,1,1,30,30,30);
            elseif isequal(obj.type,'centralSupport')
%                 obj.mesh = HexaMesh(4,4,3,32,32,24);
                obj.mesh = HexaMesh(1,1,0.75,30,30,24);
%                 obj.mesh = HexaMesh(1,1,0.75,40,40,30);
            elseif isequal(obj.type, 'torqueBeam')
%                 obj.mesh = HexaMesh(3,1,1,90,30,30);
                obj.mesh = HexaMesh(3,1,1,30,10,10);
            elseif isequal(obj.type, 'bridge')
                obj.mesh = HexaMesh(1,1,1,30,30,30);
            end
        end

        function createDesignVariable(obj)
            s.type = 'Full';
            g      = GeometricalFunction(s);
            lsFun  = g.computeLevelSetFunction(obj.mesh);
            s.fun  = lsFun;
%             s.fun.setFValues(importdata('1e-35lambda1min2gJ1cornersSupport.txt'))
%             s.fun.setFValues(importdata('1e-35lambda1min2gJ1centralSupport.txt'))
%             s.fun.setFValues(importdata('1e-35lambda1min2gJ1bridge.txt'))
%             s.fun.setFValues(importdata('1e-35lambda1min2gJ1torqueBeam.txt'))
            s.mesh = obj.mesh;
            s.type = 'LevelSet';
            s.plotting = true;
            s.isFixed.nodes = obj.createNonDesignableDomain();
            ls     = DesignVariable.create(s);
            obj.designVariable = ls;
        end

        function createFilter(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            f = Filter.create(s);
            obj.filter = f;
        end

         function createFilterCompliance(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterComp = f;
         end

        function createFilterConnectivity(obj)
            s.filterType = 'LUMP';
            s.mesh       = obj.mesh;
            s.trial      = LagrangianFunction.create(obj.mesh,1,'P1');
            f            = Filter.create(s);
            obj.filterConnect = f;
        end

        function createMaterialInterpolator(obj)
            E0   = 1e-3;
            nu0  = 1/3;
            E1   = 1;
            nu1  = 1/3;
            ndim = 3;

            matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E0,nu0);
            matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E0,nu0,ndim);

            matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E1,nu1);
            matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E1,nu1,ndim);

%             s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;
 
            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function createElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.createMaterial();
            s.dim = '3D';
            s.boundaryConditions = obj.createElasticBoundaryConditions();
            s.interpolationType = 'LINEAR';
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'CG';
%             s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                       = obj.mesh;
            s.filter                     = obj.filterComp;
            s.complainceFromConstitutive = obj.createComplianceFromConstitutive();
            s.material                   = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end
     
        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            if isequal(obj.type, 'bridge')
                s.volumeTarget = 0.4;
%             elseif isequal(obj.type, 'torqueBeam')
%                 s.volumeTarget = 0.25;
            else
                s.volumeTarget = 0.3;
            end
            v = VolumeConstraint(s);
            obj.volume = v;
        end

        function createEigenValueConstraint(obj)                           
            s.mesh              = obj.mesh;
            s.designVariable    = obj.designVariable;
            s.filter            = obj.filterConnect;
            s.filterAdjoint     = obj.filterAdjointConnect;  
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            s.targetEigenValue  = obj.lambda1min;       
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
            s.weights           = [1.0];
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
            s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            obj.constraint      = Constraint(s);
        end

%       function createDualVariable(obj)
%             s.nConstraints   = 2;
%             l                = DualVariable(s);
%             obj.dualVariable = l;
%         end  

        function createPrimalUpdater(obj)
            s.mesh = obj.mesh;
            obj.primalUpdater = SLERP(s);
        end

        function createOptimizer(obj)
            obj.gJ = 1.0;
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.type           = obj.type;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
%             s.dualVariable   = obj.dualVariable;
%             s.GIFname        = '1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.eta)+'GIF';
            s.GIFname        = '1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+'GIF';
            s.maxIter        = 1000;
            s.tolerance      = 1e-3; %3.0e-2; %
            s.constraintCase{1} = 'EQUALITY';
            s.constraintCase{2} = 'INEQUALITY';      
            s.primalUpdater  = obj.primalUpdater;
            s.etaNorm        = 0.02;  
            s.etaNormMin     = 0.02;
            s.gJFlowRatio    = 1.0;  
            s.etaMax         = 60.0; 
            s.etaMaxMin      = 1.0;  
            s.filter         = obj.filterComp;
            opt = OptimizerNullSpace(s);
            opt.solveProblem();
            obj.optimizer = opt;
            saveas(figure(1),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.type)+'design.png','png')
            saveas(figure(2),'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.type)+'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues,'1e-35lambda1min'+string(obj.lambda1min)+'gJ'+string(obj.gJ)+string(obj.type)+'.txt')
            obj.designVariable.fun.print('level_set'+string(obj.type),'Paraview')
        end

        function m = createMaterial(obj)
            x = obj.designVariable;
            f = x.obtainDomainFunction();
            f = obj.filterComp.compute(f{1},1);            
            s.type                 = 'DensityBased';
            s.density              = f;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '3D';
            s.mesh                 = obj.mesh;
            m = Material.create(s);
        end

        function bc = createElasticBoundaryConditions(obj)
            if isequal(obj.type, 'centralSupport')
                yMin    = min(obj.mesh.coord(:,2));
                xMax    = max(obj.mesh.coord(:,1));
                yMax    = max(obj.mesh.coord(:,2));
                zMax    = max(obj.mesh.coord(:,3));
                isDir   = @(coor)  (abs(coor(:,3))==0 & abs(coor(:,2))<=0.1*yMax & abs(coor(:,1))>=0.9*xMax);
                %                 isDir   = @(coor)  (abs(coor(:,3))==0 & abs(coor(:,2))<=0.25*yMax & abs(coor(:,1))>=0.75*xMax);
%                 isDir   = @(coor)  (abs(coor(:,3))==0 & abs(coor(:,2))<=0.3*yMax & abs(coor(:,1))>=0.7*xMax);
                isDirY   = @(coor)  abs(coor(:,2))==yMin;
                isDirX   = @(coor)  abs(coor(:,1))==xMax;
                isForce = @(coor)  (abs(coor(:,3))==zMax);
    
                sDir{1}.domain    = @(coor) isDir(coor);
                sDir{1}.direction = [1,2,3];
                sDir{1}.value     = 0;
    
                sDir{2}.domain    = @(coor) isDirY(coor);
                sDir{2}.direction = [2];
                sDir{2}.value     = 0;
    
                sDir{3}.domain    = @(coor) isDirX(coor);
                sDir{3}.direction = [1];
                sDir{3}.value     = 0;
    
                sPL{1}.domain    = @(coor) isForce(coor);
                sPL{1}.direction = 3;
                sPL{1}.value     = -1;
            elseif isequal(obj.type, 'torqueBeam')
                xMin    = min(obj.mesh.coord(:,1));
                xMax    = max(obj.mesh.coord(:,1));
                yMin    = min(obj.mesh.coord(:,2));
                yMax    = max(obj.mesh.coord(:,2));
                zMin    = min(obj.mesh.coord(:,3));
                zMax    = max(obj.mesh.coord(:,3));
                isDir   = @(coor)  abs(coor(:,1))==xMin;
                isForceLeft = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))==yMin  & abs(coor(:,3))>=0.45*zMax  & abs(coor(:,3))<=0.55*zMax);
                isForceDown = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,3))==zMin  & abs(coor(:,2))>=0.45*yMax  & abs(coor(:,2))<=0.55*yMax);
                isForceRight = @(coor) (abs(coor(:,1))==xMax & abs(coor(:,2))==yMax  & abs(coor(:,3))>=0.45*zMax & abs(coor(:,3))<=0.55*zMax);
                isForceUp = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,3))==zMax & abs(coor(:,2))>=0.45*yMax  & abs(coor(:,2))<=0.55*yMax);

                sDir{1}.domain    = @(coor) isDir(coor);
                sDir{1}.direction = [1,2,3];
                sDir{1}.value     = 0;
    
                sPL{1}.domain    = @(coor) isForceLeft(coor);
                sPL{1}.direction = 3;
                sPL{1}.value     = -1;
      
                sPL{2}.domain    = @(coor) isForceDown(coor);
                sPL{2}.direction = 2;
                sPL{2}.value     = 1;

                sPL{3}.domain    = @(coor) isForceRight(coor);
                sPL{3}.direction = 3;
                sPL{3}.value     = 1;
    
                sPL{4}.domain    = @(coor) isForceUp(coor);
                sPL{4}.direction = 2;
                sPL{4}.value     = -1;  
            elseif isequal(obj.type, 'cornersSupport')
                xMax    = max(obj.mesh.coord(:,1));
                yMin    = min(obj.mesh.coord(:,2));
                yMax    = max(obj.mesh.coord(:,2));
                zMin    = min(obj.mesh.coord(:,3));
                isDir2   = @(coor)  (abs(coor(:,3))==zMin & abs(coor(:,2))>=0.9*yMax & abs(coor(:,1))<=0.1*xMax);
                isForce = @(coor)  (abs(coor(:,3))==zMin & abs(coor(:,2))<=0.1*yMax & abs(coor(:,1))>=0.9*xMax);
                isDirY   = @(coor)  abs(coor(:,2))==yMin;
                isDirX   = @(coor)  abs(coor(:,1))==xMax;

                sDir{1}.domain    = @(coor) isDir2(coor);
                sDir{1}.direction = [1,2,3];
                sDir{1}.value     = 0;
                
                sDir{2}.domain    = @(coor) isDirY(coor);
                sDir{2}.direction = [2];
                sDir{2}.value     = 0;
    
                sDir{3}.domain    = @(coor) isDirX(coor);
                sDir{3}.direction = [1];
                sDir{3}.value     = 0;

                sPL{1}.domain    = @(coor) isForce(coor);
                sPL{1}.direction = 3;
                sPL{1}.value     = -1;
            elseif isequal(obj.type, 'bridge')
                yMin    = min(obj.mesh.coord(:,2));
                xMax    = max(obj.mesh.coord(:,1));
                xMin    = min(obj.mesh.coord(:,1));
                yMax    = max(obj.mesh.coord(:,2));
                zMax    = max(obj.mesh.coord(:,3));
                zMin    = min(obj.mesh.coord(:,3));
                isDir   = @(coor)  (abs(coor(:,3))==zMin & abs(coor(:,1))<=0.1*xMax);
                isDirX   = @(coor)  abs(coor(:,1))==xMax;
                isForce = @(coor)  (abs(coor(:,3))==zMax);
    
                sDir{1}.domain    = @(coor) isDir(coor);
                sDir{1}.direction = [1,2,3];
                sDir{1}.value     = 0;
    
                sDir{2}.domain    = @(coor) isDirX(coor);
                sDir{2}.direction = 1;
                sDir{2}.value     = 0;
    
                sPL{1}.domain    = @(coor) isForce(coor);
                sPL{1}.direction = 3;
                sPL{1}.value     = -1;
            end
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
    
            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
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
            zMin    = min(obj.mesh.coord(:,3));
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            zMax    = max(obj.mesh.coord(:,3));
            isDown  = @(coor) abs(coor(:,3))==zMin;
            isUp    = @(coor) abs(coor(:,3))==zMax;
            isLeft  = @(coor) abs(coor(:,1))==xMin;
            isRight = @(coor) abs(coor(:,1))==xMax;
            isFront = @(coor) abs(coor(:,2))==yMin;
            isBack = @(coor) abs(coor(:,2))== yMax;

            if isequal(obj.type, 'centralSupport') || isequal(obj.type, 'cornersSupport') 
                isDir   = @(coor)  isDown(coor) | isUp(coor) | isLeft(coor) | isBack(coor);  
            elseif isequal(obj.type, 'bridge')
                isDir   = @(coor)  (isDown(coor) | isUp(coor) | isLeft(coor) | isBack(coor) | isFront(coor));  
             else
                isDir   = @(coor)  (isDown(coor) | isUp(coor) | isLeft(coor) | isBack(coor) |isRight(coor)|isFront(coor));  
            end 

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

        function dofs = createNonDesignableDomain(obj)
            zMax    = max(obj.mesh.coord(:,3));
            xMax    = max(obj.mesh.coord(:,1));
            if isequal(obj.type, 'centralSupport')
%                 isNonDesign =  @(coor) abs(coor(:,3)) >= (zMax - 2*obj.mesh.computeMeanCellSize());  
                isNonDesign =  @(coor) abs(coor(:,3)) >= (zMax - 3*obj.mesh.computeMeanCellSize());  
            elseif isequal(obj.type, 'bridge') 
                isNonDesign =  @(coor) abs(coor(:,3)) >= (zMax - 2*obj.mesh.computeMeanCellSize());  
            elseif isequal(obj.type, 'torqueBeam')
                isNonDesign =  @(coor) ((abs(coor(:,1)) >= (xMax - 2*obj.mesh.computeMeanCellSize())) | (abs(coor(:,1)) <= (2*obj.mesh.computeMeanCellSize())));
            elseif isequal(obj.type, 'cornersSupport') 
                isNonDesign =  @(coor) abs(coor(:,3)) <= (1*obj.mesh.computeMeanCellSize());
            end
            dofs = isNonDesign(obj.mesh.coord);
        end

    end

end