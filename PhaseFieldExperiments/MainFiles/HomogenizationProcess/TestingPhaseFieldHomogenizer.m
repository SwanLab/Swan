classdef TestingPhaseFieldHomogenizer < handle

    properties (Access = private)
        E
        nu
        meshType
        meshN
        holeType
        nSteps
        damageType
        pnorm
    end

    properties (Access = private)
        baseMesh
        backgroundMesh
        test
        masterSlave
        maxParam
    end

    methods (Access = public)

        function obj = TestingPhaseFieldHomogenizer(cParams)
            obj.init(cParams);
            obj.defineMesh();
            figure()
            set(gcf, 'WindowState', 'maximized')
            drawnow
        end

        function [mat,phi,holeParams] = compute(obj)
            holeParams = obj.computeHoleParams();
            comb = table2array(combinations(holeParams{:}));
            nComb = size(comb,1);
            mat = zeros(3,3,nComb);
            phi = zeros(1,nComb);
            for i=1:nComb
                hole = comb(i,:);
                mat(:,:,i) = obj.computeHomogenization(hole);
                phi(i) = obj.computeDamageMetric(hole);
            end
            %mat = obj.assembleResults(mat);
            %phi = obj.assembleResults(phi);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.E         = cParams.E;
            obj.nu        = cParams.nu;
            obj.meshType  = cParams.meshType;
            obj.meshN     = cParams.meshN;
            obj.holeType  = cParams.holeType;
            obj.nSteps     = cParams.nSteps;
            obj.damageType = cParams.damageType;
            obj.pnorm      = cParams.pnorm;
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

            m = TriangleMesh(1,1,obj.meshN,obj.meshN);
            %m = QuadMesh(1,1,obj.meshN,obj.meshN);
            obj.backgroundMesh = m;
           % s.coord = MC.coord;
           % s.connec = MC.connec;
           % obj.masterSlave = MC.masterSlaveIndex;
            %obj.backgroundMesh = Mesh.create(s);

            obj.test = LagrangianFunction.create(obj.backgroundMesh,1,'P1');
        end

        function paramHole = computeHoleParams(obj)
            obj.maxParam = obj.computeMaxHoleParams();
            nParam = length(obj.maxParam);
            paramHole = cell(1,nParam);
            for i=1:nParam
                paramHole{i} = linspace(0,obj.maxParam(i),obj.nSteps(i));
            end
            paramHole{2} = 1e-2;
        end
        
        function maxV = computeMaxHoleParams(obj)
            switch obj.holeType
                case 'Circle'
                    maxV = 0.495;
                case 'Square'
                    maxV = 0.99;
                case 'Ellipse'
                    maxV = [0.98,0.98];
                case 'Rectangle'
                    maxV = [0.99,0.99];
                case 'SmoothHexagon'
                    maxV = 0.99;
            end
        end

        function matHomog = computeHomogenization(obj,l)
            dens = obj.createDensityLevelSet(l);
            mat  = obj.createDensityMaterial(dens);
            matHomog = obj.solveElasticMicroProblem(mat,dens);
        end

        function dens = createDensityLevelSet(obj,l)
            ls = obj.computeLevelSet(obj.backgroundMesh,l);
            sUm.backgroundMesh = obj.backgroundMesh;
            sUm.boundaryMesh   = obj.backgroundMesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);

            % ls = CharacteristicFunction.create(uMesh);
            % s.trial = obj.test;
            % s.mesh = obj.backgroundMesh;
            % f = FilterLump(s); 
            % dens = f.compute(ls,2);

            obj.baseMesh = uMesh.createFullInnerMesh('Matlab');
            
            % clf
            % obj.baseMesh.plot
            % axis off
            % title('')
            % drawnow
            % exportgraphics(gcf,'microFracture.gif','Append',true);                

            dens = LagrangianFunction.create(obj.baseMesh,1,'P1');
            fV = ones(size(dens.fValues));
            dens.setFValues(fV)
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
                    gPar.radius = l;
                case 'Square'
                    gPar.radius = l;
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
            % phiFun.plot;
            % colormap default
            lsCircle           = phiFun.fValues;
            ls = -lsCircle;
        end

        function mat = createDensityMaterial(obj,lsf)
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA = obj.createMaterial(obj.baseMesh,1e-6*obj.E,obj.nu);
            s.matB = obj.createMaterial(obj.baseMesh,obj.E,obj.nu);
            mI = MaterialInterpolator.create(s);

            x{1} = lsf;
            s.mesh                 = obj.baseMesh;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
            mat = Material.create(s);
        end

        function mat = createMaterial(obj,mesh,E,nu)
            young   = ConstantFunction.create(E,mesh);
            poisson = ConstantFunction.create(nu,mesh);
            bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(young,poisson,obj.backgroundMesh.ndim);
            shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(young,poisson);

            s.type  = 'ISOTROPIC';
            s.ptype = 'ELASTIC';
            s.mesh  = mesh;
            s.bulk  = bulk;
            s.shear = shear;
            s.ndim  = obj.backgroundMesh.ndim;
            mat     = Material.create(s);
        end

        function matHomog = solveElasticMicroProblem(obj,material,dens)
          dens.plot
          colormap (flipud(pink))

            s.mesh = obj.baseMesh;
            s.material = material;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);
            s.solverCase = 'DIRECT';
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            material.setDesignVariable({dens})
            fem.updateMaterial(material.obtainTensor())
            fem.solve();

            totVol = obj.backgroundMesh.computeVolume();
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
            %bc.updatePeriodicConditions(obj.masterSlave);
        end

        function coorRot = defineRotatedCoordinates(~,theta)
            x0 = 0.5; y0 = sqrt(1-0.5^2);
            coorRot = @(coor) feval(@(fun) fun(:,2),([cos(theta) sin(theta); sin(theta) cos(theta)]*(coor-[x0,y0])')');
        end

        function phi = computeDamageMetric(obj,l)
            max = obj.maxParam;
            switch obj.damageType
                case 'Area'
                    switch obj.holeType
                        case 'Circle'
                            phi = l^2/(max^2);
                        case 'Square'
                            phi = l^2/max^2;
                        case 'Ellipse'
                            phi = (l(1)*l(2))/(max(1)*max(2));
                        case 'Rectangle'
                            phi = (l(1)*l(2))/(max(1)*max(2));
                        case 'SmoothHexagon'
                            perimeter = 6*l;
                            apothem   = sqrt(l^2 - (l/2)^2);
                            phi = (perimeter*apothem)/(6*sqrt(3)/2);
                    end
                case 'Perimeter'
                    switch obj.holeType
                        case 'Circle'
                            phi = l/max;
                        case 'Ellipse'
                            phi = pi*(3*(l(1)+l(2))-sqrt((3*l(1)+l(2))*(l(1)+3*l(2))))/...
                                  pi*(3*(max(1)+max(2))-sqrt((3*max(1)+max(2))*(max(1)+3*max(2))));
                        case 'Square'
                            phi = l/max;
                        case 'Rectangle'
                            phi = (l(1)+l(2))/(max(1)+max(2));
                        case 'SmoothHexagon'
                            phi = 6*l;
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