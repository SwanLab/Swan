classdef TutorialLatticeNN < handle

    properties (Access = public)

        rho
        a
        b
        Chomog

    end

    properties (Access = private)

        E
        nu
        meshN
        holeType
        nSteps
        pnorm
        monitoring
        T

    end

    properties (Access = private)

        baseMesh
        masterSlave
        test

    end

    methods (Access = public)

        function obj = TutorialLatticeNN(computeData)

            if nargin < 1
                computeData = true;
            end

            obj.init();

            obj.computeParameterValues();

            if computeData
                obj.compute();
            end

        end

        function plotMicrostructures(obj,rhoValues,aValue,bValue)

            obj.computeTransformation(aValue,bValue);

            obj.defineMesh();

            for i = 1:length(rhoValues)

                rhoValue = rhoValues(i);

                dens = obj.createDensity(rhoValue);

                figure

                dens.plot

                shading interp

                colormap(flipud(pink))

                title(['\rho = ',num2str(rhoValue),...
                    ', a = ',num2str(aValue),...
                    ', b = ',num2str(bValue)])

            end

        end
        function matHomog = computeHomogenizationAt(...
                obj,rhoValue,aValue,bValue)

            obj.computeTransformation(aValue,bValue);

            obj.defineMesh();

            matHomog = obj.computeHomogenization(rhoValue);

        end
        function dataset = computeDataset(obj,rhoValues,geometryValues)

            rhoValues = rhoValues(:);

            nRho = length(rhoValues);
            nGeom = size(geometryValues,1);

            nSamples = nRho*nGeom;

            dataset = zeros(nSamples,9);

            for iGeom = 1:nGeom

                aValue = geometryValues(iGeom,1);
                bValue = geometryValues(iGeom,2);

                obj.computeTransformation(aValue,bValue);

                obj.defineMesh();

                rows = (iGeom-1)*nRho + (1:nRho);

                dataset(rows,1) = rhoValues;
                dataset(rows,2) = aValue;
                dataset(rows,3) = bValue;

                for iRho = 1:nRho

                    rhoValue = rhoValues(iRho);

                    matHomog = ...
                        obj.computeHomogenization(rhoValue);

                    dataset(rows(iRho),4:9) = ...
                        obj.tensorToVector(matHomog);

                end

            end

        end

    end

    methods (Access = private)

        function init(obj)

            obj.E = 1;
            obj.nu = 1/3;
            obj.meshN = 80;
            obj.holeType = 'Square';
            obj.pnorm = 'Inf';
            obj.nSteps = 15;
            obj.monitoring = false;

        end

        function computeParameterValues(obj)

            obj.rho = linspace(0,1,obj.nSteps);

            obj.a = 0:0.2:0.98;

            obj.b = 0.3:0.3:3;

        end
        function y = tensorToVector(~,C)

            y = [...
                C(1,1,1,1),...
                C(2,2,2,2),...
                C(1,1,2,2),...
                C(1,1,1,2),...
                C(2,2,1,2),...
                C(1,2,1,2)];

        end

        function compute(obj)

            nRho = length(obj.rho);
            nA = length(obj.a);
            nB = length(obj.b);

            mat = zeros(2,2,2,2,nRho,nA,nB);

            for iA = 1:nA

                aValue = obj.a(iA);

                for iB = 1:nB

                    bValue = obj.b(iB);

                    obj.computeTransformation(aValue,bValue);

                    obj.defineMesh();

                    for iRho = 1:nRho

                        rhoValue = obj.rho(iRho);

                        mat(:,:,:,:,iRho,iA,iB) = ...
                            obj.computeHomogenization(rhoValue);

                    end

                end

            end

            obj.Chomog = mat;

        end

        function computeTransformation(obj,aValue,bValue)

            k = 1/sqrt(1-aValue^2);

            obj.T = k*[bValue,aValue/bValue;
                       aValue*bValue,1/bValue];

        end

        function defineMesh(obj)

            s.T = obj.T;
            s.divUnit = obj.meshN;
            s.filename = '';

            MC = MeshCreator(s);

            MC.computeMeshNodes();

            s.coord = MC.coord;
            s.connec = MC.connec;

            obj.baseMesh = Mesh.create(s);

            obj.masterSlave = MC.masterSlaveIndex;

            obj.test = LagrangianFunction.create(obj.baseMesh,1,'P1');

        end

        function matHomog = computeHomogenization(obj,rho)

            dens = obj.createDensity(rho);

            C = obj.createDensityMaterial();

            matHomog = obj.solveElasticMicroProblem(C,dens);

        end

        function dens = createDensity(obj,rho)

            if rho == 0

                dens = obj.createConstantDensity(0);

            elseif rho == 1

                dens = obj.createConstantDensity(1);

            else

                l = obj.computeGeometryParameter(rho);

                dens = obj.createDensityLevelSet(l);

            end

        end

        function dens = createConstantDensity(obj,rho)

            dens = copy(obj.test);

            values = rho*ones(obj.baseMesh.nnodes,1);

            dens.setFValues(values);

        end

        function l = computeGeometryParameter(obj,rho)

            tol = 1e-4;

            scaleMin = 1e-10;
            scaleMax = 1;

            lMax = obj.createGeometryParameter(scaleMax);

            densMax = obj.createDensityLevelSet(lMax);

            rhoMax = obj.computeMaterialDensity(densMax);

            while rhoMax > rho

                scaleMax = 2*scaleMax;

                lMax = obj.createGeometryParameter(scaleMax);

                densMax = obj.createDensityLevelSet(lMax);

                rhoMax = obj.computeMaterialDensity(densMax);

            end

            for i = 1:40

                scale = 0.5*(scaleMin + scaleMax);

                l = obj.createGeometryParameter(scale);

                dens = obj.createDensityLevelSet(l);

                rhoCurrent = obj.computeMaterialDensity(dens);

                if abs(rhoCurrent-rho) < tol
                    return
                end

                if rhoCurrent > rho

                    scaleMin = scale;

                else

                    scaleMax = scale;

                end

            end

            scale = 0.5*(scaleMin + scaleMax);

            l = obj.createGeometryParameter(scale);

        end

        function l = createGeometryParameter(obj,scale)

            switch obj.holeType

                case {'Circle','Square','SmoothHexagon'}

                    l = scale;

                case {'Rectangle','Ellipse'}

                    l = [scale,scale];

            end

        end

        function rho = computeMaterialDensity(obj,dens)

            volume = obj.baseMesh.computeVolume();

            materialVolume = Integrator.compute(dens,obj.baseMesh,2);

            rho = materialVolume/volume;

        end

        function lsf = createDensityLevelSet(obj,l)

            ls = obj.computeLevelSet(obj.baseMesh,l);

            sUm.backgroundMesh = obj.baseMesh;

            sUm.boundaryMesh = obj.baseMesh.createBoundaryMesh;

            uMesh = UnfittedMesh(sUm);

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

            gPar.xCoorCenter = 0.5;
            gPar.yCoorCenter = 0.5;

            gPar.T = obj.T;

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

                    gPar.normal = [0 1;
                        sqrt(3)/2 1/2;
                        sqrt(3)/2 -1/2];

            end

            g = GeometricalFunction(gPar);

            phiFun = g.computeLevelSetFunction(mesh);

            lsGeometry = phiFun.fValues;

            ls = -lsGeometry;

        end

        function mu = computeMu(obj,E,nu)

            mu = E./(2*(1+nu));

        end

        function kappa = computeKappa(obj,E,nu,N)

            kappa = E./(N*(1-(N-1)*nu));

        end

        function C = createDensityMaterial(obj)

            N = obj.baseMesh.ndim;

            muA = obj.computeMu(1e-6*obj.E,obj.nu);

            kappaA = obj.computeKappa(1e-6*obj.E,obj.nu,N);

            muB = obj.computeMu(obj.E,obj.nu);

            kappaB = obj.computeKappa(obj.E,obj.nu,N);

            mu = @(rho) SimpAllInterpolator.computeMu(...
                muA,muB,kappaA,kappaB,rho,N);

            mu = @(rho) Expand(mu(rho),4);

            kappa = @(rho) SimpAllInterpolator.computeKappa(...
                muA,muB,kappaA,kappaB,rho,N);

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

            s.boundaryConditions = obj.createBoundaryConditions();

            s.solverCase = DirectSolver();

            s.solverType = 'REDUCED';

            s.solverMode = 'FLUC';

            fem = ElasticProblemMicro(s);

            fem.updateMaterial(C(dens));

            fem.solve();

            totVol = obj.baseMesh.computeVolume();

            matHomog = fem.Chomog/totVol;

        end

        function bc = createBoundaryConditions(obj)

            vertices = obj.baseMesh.coord(1:4,:);

            isVertex = @(coor) obj.isVertex(coor,vertices);

            sDir{1}.domain = @(coor) isVertex(coor);

            sDir{1}.direction = [1,2];

            sDir{1}.value = 0;

            dirichletFun = [];

            for i = 1:numel(sDir)

                dir = DirichletCondition(obj.baseMesh,sDir{i});

                dirichletFun = [dirichletFun,dir];

            end

            s.dirichletFun = dirichletFun;

            s.pointloadFun = [];

            s.periodicFun = 1;

            s.mesh = obj.baseMesh;

            bc = BoundaryConditions(s);

            bc.updatePeriodicConditions(obj.masterSlave);

        end

        function is = isVertex(~,coor,vertices)

            tol = 1e-12;

            is = false(size(coor,1),1);

            for i = 1:size(vertices,1)

                distance = vecnorm(coor-vertices(i,:),2,2);

                is = is | distance < tol;

            end

        end

    end

end