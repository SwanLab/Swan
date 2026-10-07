classdef TutorialHomogenizationCompetitionPlots < handle

    properties (Access = public)
        paramXi
        paramTXi
        paramA
        paramB
        Cxi
        Cab
        Csolid
    end

    properties (Access = private)
        E
        nu
        meshN
        holeType
        monitoring
        latticeVectors
        baseMesh
        masterSlave
        test
        rhoFixed
        maxParamB
        minParamB
        maxParamA
        anisotropyCoupling
        nStepsXi
        nStepsA
        nStepsB
        bFixedForAPlot
        aFixedForBPlot
    end

    methods (Access = public)

        function obj = TutorialHomogenizationCompetitionPlots()
            obj.init();
            obj.computeSolidReference();
            obj.computeXiCompetition();
            obj.computeABCompetition();
            obj.plotXiCompetition();
            obj.plotXiGeometry();
            obj.plotACompetition();
            obj.plotBCompetition();
            obj.plotABCompetition();
        end

    end

    methods (Access = private)

        function init(obj)

            obj.E = 1;
            obj.nu = 0.3;
            obj.meshN = 40;
            obj.holeType = 'RectangleAffine';
            obj.monitoring = false;

            obj.rhoFixed = 0.5;

            obj.maxParamB = 3.0;
            obj.minParamB = 1/obj.maxParamB;

            obj.maxParamA = 0.80;

            obj.anisotropyCoupling = 0.9;

            obj.nStepsXi = 31;
            obj.nStepsA  = 17;
            obj.nStepsB  = 21;

            obj.bFixedForAPlot = 1.0;
            obj.aFixedForBPlot = 0.0;

        end

        function computeSolidReference(obj)

            obj.createLattice(0,1);
            Ch = obj.computeHomogenizationFromM(0,0);
            obj.Csolid = obj.extractComponents(Ch);

        end

        function computeXiCompetition(obj)

            rho = obj.rhoFixed;
            m0 = sqrt(1-rho);

            xiMin = m0/0.98;
            xiMax = 0.98/m0;

            obj.paramXi = sort(unique([ ...
                linspace(xiMin,xiMax,obj.nStepsXi), ...
                1]));

            obj.paramTXi = ...
                (log(obj.paramXi)-log(xiMin)) ./ ...
                (log(xiMax)-log(xiMin));

            nXi = numel(obj.paramXi);
            obj.Cxi = zeros(4,nXi);

            obj.createLattice(0,1);

            fprintf('\n');
            fprintf('========================================\n');
            fprintf(' Cij(txi), rho = %.3f\n',rho);
            fprintf('========================================\n');

            for iXi = 1:nXi

                xi = obj.paramXi(iXi);

                m1 = m0*xi;
                m2 = m0/xi;

                Ch = obj.computeHomogenizationFromM(m1,m2);
                obj.Cxi(:,iXi) = obj.extractComponents(Ch);

                fprintf( ...
                    'xi = %.6f   txi = %.6f   m1 = %.6f   m2 = %.6f\n', ...
                    xi,obj.paramTXi(iXi),m1,m2);

            end

        end

        function computeABCompetition(obj)

            rho = obj.rhoFixed;

            obj.paramA = sort(unique([ ...
                linspace(0,obj.maxParamA,obj.nStepsA), ...
                0]));

            obj.paramB = sort(unique([ ...
                linspace(obj.minParamB,obj.maxParamB,obj.nStepsB), ...
                1]));

            nA = numel(obj.paramA);
            nB = numel(obj.paramB);

            obj.Cab = zeros(4,nA,nB);

            fprintf('\n');
            fprintf('========================================\n');
            fprintf(' Cij(a,b), rho = %.3f\n',rho);
            fprintf('========================================\n');

            for iA = 1:nA

                a = obj.paramA(iA);
                fprintf('\na = %.6f\n',a);

                for iB = 1:nB

                    b = obj.paramB(iB);

                    obj.createLattice(a,b);

                    [m1,m2] = obj.computeMFromBRho(b,rho);

                    Ch = obj.computeHomogenizationFromM(m1,m2);

                    obj.Cab(:,iA,iB) = obj.extractComponents(Ch);

                    fprintf( ...
                        '  b = %.6f   det(T) = %.12f\n', ...
                        b, ...
                        obj.computeDetT(a,b));

                end
            end

        end

        function createLattice(obj,a,b)

            if abs(a) >= 1
                error('Parameter a must satisfy |a| < 1.');
            end

            if b <= 0
                error('Parameter b must satisfy b > 0.');
            end

            k = 1/sqrt(1-a^2);

            B = k*[ ...
                1, a; ...
                a, 1];

            D = [ ...
                b,   0; ...
                0, 1/b];

            T = B*D;

            obj.latticeVectors = T.';
            obj.defineMesh();

        end

        function d = computeDetT(obj,a,b)

            k = 1/sqrt(1-a^2);

            B = k*[ ...
                1, a; ...
                a, 1];

            D = [ ...
                b,   0; ...
                0, 1/b];

            T = B*D;
            d = det(T);

        end

        function [m1,m2] = computeMFromBRho(obj,b,rho)

            eta = ...
                obj.anisotropyCoupling * ...
                (b - 1/b) / ...
                (obj.maxParamB - 1/obj.maxParamB);

            eta = max( ...
                -obj.anisotropyCoupling, ...
                min(obj.anisotropyCoupling,eta));

            voidFraction = 1-rho;

            m1 = ...
                voidFraction.^((1-eta)/2);

            m2 = ...
                voidFraction.^((1+eta)/2);

        end

        function Ch = computeHomogenizationFromM(obj,m1,m2)

            dens = obj.createDensityFromM1M2(m1,m2);
            mat  = obj.createDensityMaterial(dens);

            Ch = obj.solveElasticMicroProblem(mat,dens);

        end

        function lsf = createDensityFromM1M2(obj,m1,m2)

            gPar.type = obj.holeType;

            coord = obj.baseMesh.coord;

            center_x = ...
                0.5*( ...
                min(coord(:,1)) + ...
                max(coord(:,1)));

            center_y = ...
                0.5*( ...
                min(coord(:,2)) + ...
                max(coord(:,2)));

            gPar.xCoorCenter = center_x;
            gPar.yCoorCenter = center_y;

            gPar.a1 = obj.latticeVectors(1,:);
            gPar.a2 = obj.latticeVectors(2,:);

            gPar.xSide = m1;
            gPar.ySide = m2;

            g = GeometricalFunction(gPar);

            phiFun = ...
                g.computeLevelSetFunction(obj.baseMesh);

            ls = -phiFun.fValues;

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

        function c = extractComponents(obj,Ch)

            c = [ ...
                Ch(1,1,1,1); ...
                Ch(2,2,2,2); ...
                Ch(1,1,2,2); ...
                Ch(1,2,1,2)];

        end

        function defineMesh(obj)

            s.latticeVectors = obj.latticeVectors;
            s.divUnit = obj.meshN;
            s.filename = '';

            MC = MeshCreator(s);
            MC.computeMeshNodes();

            s.coord  = MC.coord;
            s.connec = MC.connec;

            obj.baseMesh = Mesh.create(s);
            obj.masterSlave = MC.masterSlaveIndex;

            obj.test = ...
                LagrangianFunction.create( ...
                obj.baseMesh, ...
                1, ...
                'P1');

        end

        function mat = createDensityMaterial(obj,lsf)

            s.interpolation = 'SIMPALL';
            s.dim = '2D';

            s.matA.bulk = ...
                IsotropicElasticMaterial. ...
                computeKappaFromYoungAndPoisson( ...
                1e-6*obj.E, ...
                obj.nu, ...
                obj.baseMesh.ndim);

            s.matA.shear = ...
                IsotropicElasticMaterial. ...
                computeMuFromYoungAndPoisson( ...
                1e-6*obj.E, ...
                obj.nu);

            s.matB.bulk = ...
                IsotropicElasticMaterial. ...
                computeKappaFromYoungAndPoisson( ...
                obj.E, ...
                obj.nu, ...
                obj.baseMesh.ndim);

            s.matB.shear = ...
                IsotropicElasticMaterial. ...
                computeMuFromYoungAndPoisson( ...
                obj.E, ...
                obj.nu);

            mI = MaterialInterpolator.create(s);

            x{1} = lsf;

            s.mesh = obj.baseMesh;
            s.type = 'DensityBased';
            s.density = x;
            s.materialInterpolator = mI;
            s.dim = '2D';

            mat = Material.create(s);

        end

        function matHomog = solveElasticMicroProblem(obj,material,dens)

            if obj.monitoring
                close all
                dens.plot
                shading interp
                colormap(flipud(pink))
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

            material.setDesignVariable({dens});
            fem.updateMaterial(material.obtainTensor());
            fem.solve();

            totVol = obj.baseMesh.computeVolume();
            matHomog = fem.Chomog/totVol;

        end

        function bc = createBoundaryConditions(obj,mesh)

            coord = mesh.coord;
            cornerCoord = coord(1:4,:);
            tol2 = 1e-20;

            isCorner = @(coor) ...
                (sum((coor-cornerCoord(1,:)).^2,2) < tol2) | ...
                (sum((coor-cornerCoord(2,:)).^2,2) < tol2) | ...
                (sum((coor-cornerCoord(3,:)).^2,2) < tol2) | ...
                (sum((coor-cornerCoord(4,:)).^2,2) < tol2);

            sDir{1}.domain = @(coor) isCorner(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value = 0;

            dirichletFun = [];

            for i = 1:numel(sDir)
                dirichletFun = [ ...
                    dirichletFun, ...
                    DirichletCondition(mesh,sDir{i})];
            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = 1;
            s.mesh = mesh;

            bc = BoundaryConditions(s);
            bc.updatePeriodicConditions(obj.masterSlave);

        end

        function plotXiCompetition(obj)

            names = { ...
                'C_{11}', ...
                'C_{22}', ...
                'C_{12}', ...
                'C_{33}'};

            [~,iRef] = min(abs(obj.paramXi-1));

            Cref = obj.Cxi(:,iRef);

            Cnorm = obj.Cxi ./ Cref;

            figure( ...
                'Color','w', ...
                'Position',[100 100 1200 850]);

            t = tiledlayout( ...
                2,2, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            for i = 1:4

                nexttile

                plot( ...
                    obj.paramXi, ...
                    Cnorm(i,:), ...
                    'LineWidth',2);

                hold on

                xline(1,':');
                yline(1,':');

                xlabel('\xi');

                ylabel('C_{IJ}/C_{IJ}(\xi=1)');

                title(names{i});

                xlim([min(obj.paramXi) max(obj.paramXi)]);

                grid on
                box on

            end

            title( ...
                t, ...
                sprintf( ...
                'Normalized C_{IJ}(\\xi), \\rho = %.2f', ...
                obj.rhoFixed));

        end

        function plotXiGeometry(obj)

            rho = obj.rhoFixed;

            m0 = sqrt(1-rho);

            m1 = m0*obj.paramXi;
            m2 = m0./obj.paramXi;

            figure( ...
                'Color','w', ...
                'Position',[150 150 850 650]);

            plot( ...
                obj.paramXi, ...
                m1, ...
                'LineWidth',2);

            hold on

            plot( ...
                obj.paramXi, ...
                m2, ...
                '--', ...
                'LineWidth',2);

            xline(1,':');

            xlabel('\xi');
            ylabel('m');

            legend( ...
                'm_1', ...
                'm_2', ...
                'Location','best');

            title( ...
                sprintf( ...
                'm_1 and m_2 vs \\xi, \\rho = %.2f', ...
                rho));

            xlim([min(obj.paramXi) max(obj.paramXi)]);
            ylim([0 1]);

            grid on
            box on

        end

        function plotACompetition(obj)

            names = { ...
                'C_{11}', ...
                'C_{22}', ...
                'C_{12}', ...
                'C_{33}'};

            [~,iA0] = min(abs(obj.paramA));
            [~,iB1] = min(abs(obj.paramB-1));

            Cref = obj.Cab(:,iA0,iB1);

            data = squeeze( ...
                obj.Cab(:,:,iB1));

            dataNorm = ...
                data ./ Cref;

            figure( ...
                'Color','w', ...
                'Position',[120 120 1200 850]);

            t = tiledlayout( ...
                2,2, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            for i = 1:4

                nexttile

                plot( ...
                    obj.paramA, ...
                    dataNorm(i,:), ...
                    'LineWidth',2);

                hold on

                xline(0,':');
                yline(1,':');

                xlabel('a');

                ylabel( ...
                    'C_{IJ}/C_{IJ}(a=0,b=1)');

                title(names{i});

                grid on
                box on

            end

            title( ...
                t, ...
                sprintf( ...
                ['Normalized C_{IJ}(a), ' ...
                '\\rho = %.2f, b = 1'], ...
                obj.rhoFixed));

        end

        function plotBCompetition(obj)

            names = { ...
                'C_{11}', ...
                'C_{22}', ...
                'C_{12}', ...
                'C_{33}'};

            [~,iA0] = min(abs(obj.paramA));
            [~,iB1] = min(abs(obj.paramB-1));

            Cref = obj.Cab(:,iA0,iB1);

            data = squeeze( ...
                obj.Cab(:,iA0,:));

            dataNorm = ...
                data ./ Cref;

            figure( ...
                'Color','w', ...
                'Position',[140 140 1200 850]);

            t = tiledlayout( ...
                2,2, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            for i = 1:4

                nexttile

                plot( ...
                    obj.paramB, ...
                    dataNorm(i,:), ...
                    'LineWidth',2);

                hold on

                xline(1,':');
                yline(1,':');

                xlabel('b');

                ylabel( ...
                    'C_{IJ}/C_{IJ}(a=0,b=1)');

                title(names{i});

                grid on
                box on

            end

            title( ...
                t, ...
                sprintf( ...
                ['Normalized C_{IJ}(b), ' ...
                '\\rho = %.2f, a = 0'], ...
                obj.rhoFixed));

        end

        function plotABCompetition(obj)

            names = { ...
                'C_{11}', ...
                'C_{22}', ...
                'C_{12}', ...
                'C_{33}'};

            [~,iA0] = min(abs(obj.paramA));
            [~,iB1] = min(abs(obj.paramB-1));

            Cref = obj.Cab(:,iA0,iB1);

            [B,A] = ...
                meshgrid( ...
                obj.paramB, ...
                obj.paramA);

            figure( ...
                'Color','w', ...
                'Position',[160 100 1250 900]);

            t = tiledlayout( ...
                2,2, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            for i = 1:4

                nexttile

                data = squeeze( ...
                    obj.Cab(i,:,:));

                dataNorm = ...
                    data / Cref(i);

                surf( ...
                    A, ...
                    B, ...
                    dataNorm, ...
                    'EdgeColor','none');

                hold on

                plot3( ...
                    0, ...
                    1, ...
                    1, ...
                    'ko', ...
                    'MarkerFaceColor','k', ...
                    'MarkerSize',7);

                xlabel('a');
                ylabel('b');

                zlabel( ...
                    'C_{IJ}/C_{IJ}(a=0,b=1)');

                title(names{i});

                grid on
                box on
                colorbar

                view(45,30);

            end

            title( ...
                t, ...
                sprintf( ...
                ['Normalized C_{IJ}(a,b), ' ...
                '\\rho = %.2f'], ...
                obj.rhoFixed));

        end

    end

end