classdef TutorialXXPhaseFieldHomogenization < handle

    properties (Access = public)
        paramHole
        Chomog
        damage
    end

    properties (Access = private)
        E
        nu
        meshType
        meshN
        holeType
        nSteps
        damageType
        pnorm
        monitoring
    end

    properties (Access = private)
        baseMesh
        masterSlave
        test
        maxParam
    end

    methods (Access = public)
        
        function obj = TutorialXXPhaseFieldHomogenization()
            obj.init();
            obj.defineMesh();
            obj.computeHoleParams();
            obj.compute();
        end
      
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.E          = 210;
            obj.nu         = 0.3;
            obj.meshType   = 'Square';
            obj.meshN      = 100;

            obj.holeType   = 'Square';
            obj.pnorm      = 'Inf';
            obj.damageType = 'Area';
            obj.nSteps     = 5;

            obj.monitoring = true;
        end

        function defineMesh(obj)
            switch obj.meshType
                case 'Square'
                    s.c = [1,1];
                    s.theta = [0,90];
                    s.divUnit = obj.meshN;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
                case 'Hexagon'
                    s.c = [1,1,1];
                    s.theta = [0,60,120];
                    s.divUnit = obj.meshN;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
            end
            s.coord  = MC.coord;
            s.connec = MC.connec;
            obj.baseMesh = Mesh.create(s);
            obj.masterSlave = MC.masterSlaveIndex;
            obj.test = LagrangianFunction.create(obj.baseMesh,1,'P1');
        end

        function computeHoleParams(obj)
            obj.maxParam = 0.979*ones(size(obj.nSteps));
            nParam = length(obj.maxParam);
            obj.paramHole = cell(1,nParam);
            for i=1:nParam
                obj.paramHole{i} = linspace(1e-5,obj.maxParam(i),obj.nSteps(i));
            end
        end

        function compute(obj)
            comb = table2array(combinations(obj.paramHole{:}));
            nComb = size(comb,1);
            mat = zeros(2,2,2,2,nComb);
            phi = zeros(1,nComb);
            for i=1:nComb
                hole = comb(i,:);
                if i==1
                    hole = 1e-10*ones(size(hole));
                end
                mat(:,:,:,:,i) = obj.computeHomogenization(hole);
                phi(i)     = obj.computeDamageMetric(hole);
            end
            obj.Chomog = obj.assembleResults(mat);
            obj.damage = obj.assembleResults(phi);
        end
        
        function matHomog = computeHomogenization(obj,l)
            dens = obj.createDensityLevelSet(l);
            mat  = obj.createDensityMaterial(dens);
            matHomog = obj.solveElasticMicroProblem(mat,dens);
        end

        function lsf = createDensityLevelSet(obj,l)
            ls = obj.computeLevelSet(obj.baseMesh,l);
            sUm.backgroundMesh = obj.baseMesh;
            sUm.boundaryMesh   = obj.baseMesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);

            ls = CharacteristicFunction.create(uMesh);
            s.trial = obj.test;
            s.mesh = obj.baseMesh;
            f = FilterLump(s); 
            lsf = f.compute(ls,2);
        end

        function ls = computeLevelSet(obj,mesh,l)
            gPar.type = obj.holeType;
            gPar.pnorm = obj.pnorm;
            switch obj.meshType
                case 'Square'
                    gPar.xCoorCenter = 0.5;
                    gPar.yCoorCenter = 0.5;
                case 'Hexagon'
                    gPar.xCoorCenter = 0.5;
                    gPar.yCoorCenter = sqrt(1-0.5^2);
            end
            switch obj.holeType
                case 'Circle'
                    gPar.radius = l/2;
                case 'Square'
                    gPar.length = l;
                case 'Rectangle'
                    gPar.xSide  = l(1);
                    gPar.ySide  = l(2);
                case 'Ellipse'
                    gPar.type = "SmoothRectangle";
                    gPar.xSide  = l(1);
                    gPar.ySide  = l(2);
                    gPar.pnorm  = 2;  
                case 'SmoothHexagon'
                    gPar.radius = l;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
            end
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            ls = -lsCircle;
        end

        function mat = createDensityMaterial(obj,lsf)
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1e-6*obj.E,obj.nu,obj.baseMesh.ndim);
            s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1e-6*obj.E,obj.nu);
            s.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(obj.E,obj.nu,obj.baseMesh.ndim);
            s.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(obj.E,obj.nu);
            mI = MaterialInterpolator.create(s);

            x{1} = lsf;
            s.mesh                 = obj.baseMesh;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
            mat = Material.create(s);
        end

        function matHomog = solveElasticMicroProblem(obj,material,dens)
            if obj.monitoring == true
                close all
                dens.plot
                shading interp
                colormap (flipud(pink))
                drawnow
            end

            s.mesh = obj.baseMesh;
            s.material = material;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            material.setDesignVariable({dens})
            fem.updateMaterial(material.obtainTensor())
            fem.solve();

            totVol = obj.baseMesh.computeVolume();
            matHomog = fem.Chomog/totVol;
        end

        function bc = createBoundaryConditions(obj,mesh)
            switch obj.meshType
                case 'Square'
                    isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
                    isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
                    isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1))) < 1e-12);
                    isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12);
                    isVertex = @(coor) (isTop(coor) & isLeft(coor))    |...
                                       (isTop(coor) & isRight(coor))   |...
                                       (isBottom(coor) & isLeft(coor)) |...
                                       (isBottom(coor) & isRight(coor));
                    sDir{1}.domain    = @(coor) isVertex(coor);
                    sDir{1}.direction = [1,2];
                    sDir{1}.value     = 0;
                case 'Hexagon'
                    isBottom      = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
                    isTop         = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
                    
                    coorRotY = obj.defineRotatedCoordinates(pi/3);
                    isRightBottom = @(coor) (abs(coorRotY(coor) - min(coorRotY(coor))) < 1e-12);
                    isLeftTop     = @(coor) (abs(coorRotY(coor) - max(coorRotY(coor))) < 1e-12);
                    coorRotY = obj.defineRotatedCoordinates(-pi/3);
                    isLeftBottom  = @(coor) (abs(coorRotY(coor) - min(coorRotY(coor))) < 1e-12);
                    isRightTop    = @(coor) (abs(coorRotY(coor) - max(coorRotY(coor))) < 1e-12);
                    isVertex = @(coor) (isBottom(coor) & isRightBottom(coor))  |...
                                       (isRightBottom(coor) & isRightTop(coor))|...
                                       (isRightTop(coor) & isTop(coor))        |...
                                       (isTop(coor) & isLeftTop(coor))         |...
                                       (isLeftTop(coor) & isLeftBottom(coor))  |...
                                       (isLeftBottom(coor) & isBottom(coor))   ;
                    sDir{1}.domain    = @(coor) isVertex(coor);
                    sDir{1}.direction = [1,2];
                    sDir{1}.value     = 0;
            end

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = 1; %Set to not be empty
            s.mesh = mesh;
            bc = BoundaryConditions(s);
            bc.updatePeriodicConditions(obj.masterSlave);
        end

        function coorRot = defineRotatedCoordinates(~,theta)
            x0 = 0.5; y0 = sqrt(1-0.5^2);
            coorRot = @(coor) feval(@(fun) fun(:,2),([cos(theta) sin(theta); sin(theta) cos(theta)]*(coor-[x0,y0])')');
        end

        function phi = computeDamageMetric(obj,l)
            switch obj.damageType
                case 'Area'
                    switch obj.holeType
                        case {'Circle','Square'}
                            phi = l^2;
                        case {'Ellipse','Rectangle'}
                            phi = l(1)*l(2);
                        case {'SmoothHexagon','Hexagon'}
                            perimeter = 6*l;
                            apothem   = sqrt(l^2 - (l/2)^2);
                            phi = (perimeter*apothem)/(6*sqrt(3)/2);
                    end
                case 'Perimeter'
                    switch obj.holeType
                        case {'Circle','Square','SmoothHexagon','Hexagon'}
                            phi = l;
                        case 'Ellipse'
                            phi = pi*(3*(l(1)+l(2))-sqrt((3*l(1)+l(2))*(l(1)+3*l(2))))/...
                                  pi*(3*(2)-sqrt((3+1)*(1+3)));
                        case 'Rectangle'
                            phi = (l(1)+l(2))/2;
                    end
            end
        end
        
        function [mat] = assembleResults(obj,vec)
            sizeRes = size(vec);
            mat = zeros([sizeRes(1:end-1),obj.nSteps]);
            nStepsLastParam = obj.nSteps(end);
            nCombs = sizeRes(end);
            idxVec = repmat({':'}, 1, ndims(vec));
            idxMat = repmat({':'}, 1, ndims(mat));
            for i=1:nStepsLastParam
                idxVec{end} = i:nStepsLastParam:nCombs;
                idxMat{end} = i;
                mat(idxMat{:}) = vec(idxVec{:});
            end
        end
        
    end
    
end