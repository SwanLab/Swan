classdef TopOptDensity3DConnectivity < handle

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
        perimeter
        gJ
        lambda1min
        dofsNonDesign
        type
    end

    methods (Access = public)

        function obj = TopOptDensity3DConnectivity()    
            for type = ['bridge',  'cornersSupport', 'centralSupport', 'torqueBeam']
                obj.type = type;
                obj.init()
                obj.createMesh();
                obj.createDesignVariable();
                obj.createFilterCompliance();
                obj.createFilterConnectivity();
                obj.createMaterialInterpolator();
                obj.createElasticProblem();
                obj.createComplianceFromConstitutive();
                obj.createNonDesignableDomain();
                obj.createCompliance();
                obj.createPerimeter();
                obj.createEigenValueConstraint();                             
                obj.createVolumeConstraint();
                obj.createCost();
                obj.createConstraint();
                obj.createDualVariable();
                obj.createOptimizer();
            end
         end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            if isequal(obj.type, 'cornersSupport') 
                obj.mesh = HexaMesh(2,2,2,30,30,30);
            elseif isequal(obj.type,'centralSupport')
                obj.mesh = HexaMesh(4,4,3,40,40,30);
            elseif isequal(obj.type, 'torqueBeam')
                obj.mesh = HexaMesh(3,1,1,75,25,25);
            elseif isequal(obj.type, 'bridge')
                obj.mesh = HexaMesh(1,1,1,20,20,20);
            end
        end
    
        function createDesignVariable(obj)
            s.fHandle = @(x) 1.0*ones(size(x(1,:,:)));
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            sD.fun      = aFun.project('P1');
            sD.mesh     = obj.mesh;
            sD.type     = 'Density';
            sD.plotting = true;
            dens        = DesignVariable.create(sD);
            obj.designVariable = dens;
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
            s.dim            = '3D';
            s.matA = matA;
            s.matB = matB;

            m = MaterialInterpolator.create(s);
            obj.materialInterpolator = m;
        end

        function m = createMaterial(obj)
            f = obj.designVariable.fun;           
            s.type                 = 'DensityBased';
            s.density              = f;
            s.mesh                 = obj.mesh;
            s.materialInterpolator = obj.materialInterpolator;
            s.dim                  = '3D';
            m = Material.create(s);
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
            fem = ElasticProblem(s);
            obj.physicalProblem = fem;
        end

        function c = createComplianceFromConstitutive(obj)
            s.mesh         = obj.mesh;
            s.stateProblem = obj.physicalProblem;
            c = ComplianceFromConstitutiveTensor(s);
        end

        function createCompliance(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterComp;
            if ~isempty(obj.filterAdjointComp)
                s.filterAdjoint               = obj.filterAdjointComp;
            end
            s.complainceFromConstitutive  = obj.createComplianceFromConstitutive();
            s.material                    = obj.createMaterial();
            c = ComplianceFunctional(s);
            obj.compliance = c;
        end

        function createPerimeter(obj)
            s.mesh                        = obj.mesh;
            s.filter                      = obj.filterComp;
            p = SimplePerimeterFunctional(s);
            obj.perimeter = p;
        end

        function createVolumeConstraint(obj)
            s.mesh   = obj.mesh;
            s.filter = obj.filterComp;
            if ~isempty(obj.filterAdjointComp)
                s.filterAdjoint               = obj.filterAdjointComp;
            end
            s.gradientTest = LagrangianFunction.create(obj.mesh,1,'P1');
            if isequal(obj.type, 'bridge')
                s.volumeTarget = 0.4;
            elseif isequal(obj.type, 'torqueBeam')
                s.volumeTarget = 0.4;
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
            s.targetEigenValue  = 2.0;      
            s.boundaryConditions = obj.createEigenvalueBoundaryConditions();
            obj.minimumEigenValue = StiffnessEigenModesConstraint(s);
        end

        function createConstraint(obj)
            s.shapeFunctions{1} = obj.volume;
%             s.shapeFunctions{2} = obj.minimumEigenValue; 
            s.Msmooth           = obj.createMassMatrix();
            s.dofsNonDesign     = obj.dofsNonDesign;
            obj.constraint      = Constraint(s);
        end

        function createCost(obj)
            s.shapeFunctions{1} = obj.compliance;
%             s.shapeFunctions{2} = obj.perimeter;
            s.weights           = [1.0];
            s.Msmooth           = obj.createMassMatrix();
            s.dofsNonDesign     = obj.dofsNonDesign;
            obj.cost            = Cost(s);
        end

        function M = createMassMatrix(obj)
            n = obj.mesh.nnodes;
            h = obj.mesh.computeMinCellSize();
            M = h^2*sparse(1:n,1:n,ones(1,n),n,n);
        end

        function createDualVariable(obj)
            s.nConstraints   = 1;                                        
            l                = DualVariable(s);
            obj.dualVariable = l;
        end

        function createOptimizer(obj)
            s.monitoring     = true;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.maxIter        = 1000;
            s.tolerance      = 1e-8;
            s.constraintCase{1} = 'EQUALITY';
%             s.constraintCase{2} = 'INEQUALITY';                             
            s.ub             = 1;
            s.lb             = 0;
            s.dofsNonDesign  = obj.dofsNonDesign;
            opt              = OptimizerMMA(s);
            opt.solveProblem();
            obj.optimizer = opt;
            saveas(figure(1),'DEN1e-3lambda1min'+string(obj.lambda1min)+'design.png','png')
            saveas(figure(2),'DEN1e-3lambda1min'+string(obj.lambda1min)+'graficos.png','png')
            writematrix(obj.designVariable.fun.fValues,'DEN1e-3lambda1min'+string(obj.lambda1min)+'.txt')
            obj.designVariable.fun.print('density','Paraview')
        end

        function bc = createElasticBoundaryConditions(obj)
            if isequal(obj.type, 'centralSupport')
                yMin    = min(obj.mesh.coord(:,2));
                xMax    = max(obj.mesh.coord(:,1));
                yMax    = max(obj.mesh.coord(:,2));
                zMax    = max(obj.mesh.coord(:,3));
                isDir   = @(coor)  (abs(coor(:,3))==0 & abs(coor(:,2))<=0.1*yMax & abs(coor(:,1))>=0.9*xMax);
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
                sDir{2}.direction = [1];
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
                isDir   = @(coor)  (isDown(coor) | isUp(coor) | isLeft(coor) | isBack(coor) |isRight(coor)| isFront(coor));  
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

        function createNonDesignableDomain(obj)
            zMax    = max(obj.mesh.coord(:,3));
            xMax    = max(obj.mesh.coord(:,1));
            if isequal(obj.type, 'centralSupport')
                isNonDesign =  @(coor) abs(coor(:,3)) >= (zMax - 3*obj.mesh.computeMeanCellSize());  
            elseif isequal(obj.type, 'bridge') 
                isNonDesign =  @(coor) abs(coor(:,3)) >= (zMax - 2*obj.mesh.computeMeanCellSize());  
            elseif isequal(obj.type, 'torqueBeam')
                isNonDesign =  @(coor) ((abs(coor(:,1)) >= (xMax - 2*obj.mesh.computeMeanCellSize())) | (abs(coor(:,1)) <= (2*obj.mesh.computeMeanCellSize())));
            elseif isequal(obj.type, 'cornersSupport') 
                isNonDesign =  @(coor) abs(coor(:,3)) <= (2*obj.mesh.computeMeanCellSize());
            end
%             obj.dofsNonDesign = isNonDesign(obj.mesh.coord);
        end
     end
end
