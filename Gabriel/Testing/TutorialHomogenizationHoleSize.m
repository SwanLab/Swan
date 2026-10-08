classdef TutorialHomogenizationHoleSize < handle

    properties (Access = public)

        rho
        Chomog

    end

    properties (Access = private)

        E
        nu
        meshType
        meshN
        holeType
        nSteps
        pnorm
        monitoring

    end

    properties (Access = private)

        baseMesh
        masterSlave
        test

    end

    methods (Access = public)

        function obj = TutorialHomogenizationHoleSize()

            obj.init();
            obj.defineMesh();
            obj.computeDensityParams();
            obj.compute();

        end

    end

    methods (Access = private)

        function init(obj)

            obj.E          = 1;
            obj.nu         = 1/3;
            obj.meshType   = 'Square';
            obj.meshN      = 80;
            obj.holeType   = 'Square';
            obj.pnorm      = 'Inf';
            obj.nSteps     = 80;
            obj.monitoring = false;

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

        function computeDensityParams(obj)

            obj.rho = linspace(1e-5,0.979,obj.nSteps);

        end

        function compute(obj)

            nRho = length(obj.rho);

            mat = zeros(2,2,2,2,nRho);

            for i = 1:nRho

                rho = obj.rho(i);

                mat(:,:,:,:,i) = obj.computeHomogenization(rho);

            end

            obj.Chomog = mat;

        end

        function matHomog = computeHomogenization(obj,rho)

            dens = obj.createDensityLevelSet(rho);

            C = obj.createDensityMaterial();

            matHomog = obj.solveElasticMicroProblem(C,dens);

        end

        

        function lsf = createDensityLevelSet(obj,l)

            ls = obj.computeLevelSet(obj.baseMesh,l);

            sUm.backgroundMesh = obj.baseMesh;
            sUm.boundaryMesh   = obj.baseMesh.createBoundaryMesh;

            uMesh = UnfittedMesh(sUm);

            uMesh.compute(ls);

            ls = CharacteristicFunction.create(uMesh);

            s.trial = obj.test;
            s.mesh  = obj.baseMesh;

            f = FilterLump(s);

            lsf = f.compute(ls,2);

        end

        function ls = computeLevelSet(obj,mesh,l)

            gPar.type  = obj.holeType;
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

                    gPar.xSide = l(1);
                    gPar.ySide = l(2);

                case 'Ellipse'

                    gPar.type = "SmoothRectangle";
                    gPar.xSide = l(1);
                    gPar.ySide = l(2);
                    gPar.pnorm = 2;

                case 'SmoothHexagon'

                    gPar.radius = l;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];

            end

            g = GeometricalFunction(gPar);

            phiFun = g.computeLevelSetFunction(mesh);

            lsCircle = phiFun.fValues;

            ls = -lsCircle;

        end

        function mu = computeMu(obj,E,nu)

            mu = E./(2*(1+nu));

        end

        function kappa = computeKappa(obj,E,nu,N)

            kappa = E./(N*(1-(N-1)*nu));

        end

        function C = createDensityMaterial(obj)

            N = obj.baseMesh.ndim;

            muA    = obj.computeMu(1e-6*obj.E,obj.nu);
            kappaA = obj.computeKappa(1e-6*obj.E,obj.nu,N);

            muB    = obj.computeMu(obj.E,obj.nu);
            kappaB = obj.computeKappa(obj.E,obj.nu,N);

            mu = @(rho) SimpAllInterpolator.computeMu(muA,muB,kappaA,kappaB,rho,N);

            mu = @(rho) Expand(mu(rho),4);

            kappa = @(rho) SimpAllInterpolator.computeKappa(muA,muB,kappaA,kappaB,rho,N);

            kappa = @(rho) Expand(kappa(rho),4);

            lambda = @(rho) kappa(rho) - (2/N)*mu(rho);

            I = ConstantFunction.create(eye4D(N),obj.baseMesh);

            IxI = ConstantFunction.create(kronEye(N),obj.baseMesh);

            C = @(rho) 2*mu(rho).*I + lambda(rho).*IxI;

        end

        function matHomog = solveElasticMicroProblem(obj,C,dens)

            if obj.monitoring == true

                close all

                dens.plot

                shading interp

                colormap(flipud(pink))

                drawnow

            end

            s.mesh = obj.baseMesh;

            s.material = [];

            s.scale = 'MICRO';

            s.dim = '2D';

            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);

            s.solverCase = DirectSolver();

            s.solverType = 'REDUCED';

            s.solverMode = 'FLUC';

            fem = ElasticProblemMicro(s);

            fem.updateMaterial(C(dens));

            fem.solve();

            totVol = obj.baseMesh.computeVolume();

            matHomog = fem.Chomog/totVol;

        end

        function bc = createBoundaryConditions(obj,mesh)

            switch obj.meshType

                case 'Square'

                    isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);

                    isTop = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);

                    isRight = @(coor) (abs(coor(:,1) - max(coor(:,1))) < 1e-12);

                    isLeft = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12);

                    isVertex = @(coor) (isTop(coor) & isLeft(coor)) |...
                                       (isTop(coor) & isRight(coor)) |...
                                       (isBottom(coor) & isLeft(coor)) |...
                                       (isBottom(coor) & isRight(coor));

                    sDir{1}.domain = @(coor) isVertex(coor);

                    sDir{1}.direction = [1,2];

                    sDir{1}.value = 0;

                case 'Hexagon'

                    isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);

                    isTop = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);

                    coorRotY = obj.defineRotatedCoordinates(pi/3);

                    isRightBottom = @(coor) (abs(coorRotY(coor) - min(coorRotY(coor))) < 1e-12);

                    isLeftTop = @(coor) (abs(coorRotY(coor) - max(coorRotY(coor))) < 1e-12);

                    coorRotY = obj.defineRotatedCoordinates(-pi/3);

                    isLeftBottom = @(coor) (abs(coorRotY(coor) - min(coorRotY(coor))) < 1e-12);

                    isRightTop = @(coor) (abs(coorRotY(coor) - max(coorRotY(coor))) < 1e-12);

                    isVertex = @(coor) (isBottom(coor) & isRightBottom(coor)) |...
                                       (isRightBottom(coor) & isRightTop(coor)) |...
                                       (isRightTop(coor) & isTop(coor)) |...
                                       (isTop(coor) & isLeftTop(coor)) |...
                                       (isLeftTop(coor) & isLeftBottom(coor)) |...
                                       (isLeftBottom(coor) & isBottom(coor));

                    sDir{1}.domain = @(coor) isVertex(coor);

                    sDir{1}.direction = [1,2];

                    sDir{1}.value = 0;

            end

            dirichletFun = [];

            for i = 1:numel(sDir)

                dir = DirichletCondition(mesh,sDir{i});

                dirichletFun = [dirichletFun,dir];

            end

            s.dirichletFun = dirichletFun;

            s.pointloadFun = [];

            s.periodicFun = 1;

            s.mesh = mesh;

            bc = BoundaryConditions(s);

            bc.updatePeriodicConditions(obj.masterSlave);

        end

        function coorRot = defineRotatedCoordinates(~,theta)

            x0 = 0.5;

            y0 = sqrt(1-0.5^2);

            coorRot = @(coor) feval(@(fun) fun(:,2),...
                ([cos(theta) sin(theta); sin(theta) cos(theta)]*(coor-[x0,y0])')');

        end

    end

end