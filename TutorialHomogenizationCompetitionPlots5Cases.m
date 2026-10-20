classdef TutorialHomogenizationCompetitionPlots5Cases < handle

    properties (Access = public)
        rhoList
        paramXi
        paramA
        paramB
        paramBeta
        Cxi
        Ca
        Cb
        Cab
        Cbeta
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
        maxParamB
        minParamB
        maxParamA
        nStepsXi
        nStepsA
        nStepsB
        nStepsBeta
    end

    methods (Access = public)

        function obj = TutorialHomogenizationCompetitionPlots5Cases()
            obj.init();
            obj.computeSolidReference();
            obj.computeXiCompetition();
            obj.computeACompetition();
            obj.computeBCompetition();
            obj.computeABCompetition();
            obj.computeBetaCompetition();
            obj.plotXiCompetition();
            obj.plotACompetition();
            obj.plotBCompetition();
            obj.plotABCompetition();
            obj.plotBetaCompetition();
        end

    end

    methods (Access = private)

        function init(obj)

            obj.E = 1;
            obj.nu = 0.3;
            obj.meshN = 40;
            obj.holeType = 'RectangleAffine';
            obj.monitoring = false;

            obj.rhoList = [0.2 0.5 0.8];

            obj.maxParamB = 3.0;
            obj.minParamB = 1/obj.maxParamB;

            obj.maxParamA = 0.80;

            obj.nStepsXi   = 31;
            obj.nStepsA    = 17;
            obj.nStepsB    = 21;
            obj.nStepsBeta = 21;

        end

        function computeSolidReference(obj)

            obj.createLattice(0,1);

            K = ...
                IsotropicElasticMaterial. ...
                computeKappaFromYoungAndPoisson( ...
                obj.E, ...
                obj.nu, ...
                obj.baseMesh.ndim);

            mu = ...
                IsotropicElasticMaterial. ...
                computeMuFromYoungAndPoisson( ...
                obj.E, ...
                obj.nu);

            obj.Csolid = [ ...
                K + mu; ...
                K + mu; ...
                K - mu; ...
                mu];

        end

        function computeXiCompetition(obj)

            nRho = numel(obj.rhoList);
            obj.paramXi = cell(1,nRho);
            obj.Cxi = cell(1,nRho);

            obj.createLattice(0,1);

            for iRho = 1:nRho

                rho = obj.rhoList(iRho);
                m0  = sqrt(1-rho);

                xiMin = m0/0.98;
                xiMax = 0.98/m0;

                xiVec = sort(unique([ ...
                    linspace(xiMin,xiMax,obj.nStepsXi), ...
                    1]));

                nXi = numel(xiVec);
                Ctmp = zeros(4,nXi);

                for iXi = 1:nXi

                    xi = xiVec(iXi);

                    m1 = m0*xi;
                    m2 = m0/xi;

                    Ch = obj.computeHomogenizationFromM(m1,m2);
                    Ctmp(:,iXi) = obj.extractComponents(Ch);

                end

                obj.paramXi{iRho} = xiVec;
                obj.Cxi{iRho} = Ctmp;

            end

        end

        function computeACompetition(obj)

            obj.paramA = sort(unique([ ...
                linspace(0,obj.maxParamA,obj.nStepsA), ...
                0]));

            nA   = numel(obj.paramA);
            nRho = numel(obj.rhoList);

            obj.Ca = zeros(4,nA,nRho);

            for iRho = 1:nRho

                rho = obj.rhoList(iRho);
                m0  = sqrt(1-rho);
                m1  = m0;
                m2  = m0;

                for iA = 1:nA

                    a = obj.paramA(iA);

                    obj.createLattice(a,1);

                    Ch = obj.computeHomogenizationFromM(m1,m2);
                    obj.Ca(:,iA,iRho) = obj.extractComponents(Ch);

                end
            end

        end

        function computeBCompetition(obj)

            obj.paramB = sort(unique([ ...
                linspace(obj.minParamB,obj.maxParamB,obj.nStepsB), ...
                1]));

            nB   = numel(obj.paramB);
            nRho = numel(obj.rhoList);

            obj.Cb = zeros(4,nB,nRho);

            for iRho = 1:nRho

                rho = obj.rhoList(iRho);
                m0  = sqrt(1-rho);
                m1  = m0;
                m2  = m0;

                for iB = 1:nB

                    b = obj.paramB(iB);

                    obj.createLattice(0,b);

                    Ch = obj.computeHomogenizationFromM(m1,m2);
                    obj.Cb(:,iB,iRho) = obj.extractComponents(Ch);

                end
            end

        end

        function computeABCompetition(obj)

            obj.paramA = sort(unique([ ...
                linspace(0,obj.maxParamA,obj.nStepsA), ...
                0]));

            obj.paramB = sort(unique([ ...
                linspace(obj.minParamB,obj.maxParamB,obj.nStepsB), ...
                1]));

            nA   = numel(obj.paramA);
            nB   = numel(obj.paramB);
            nRho = numel(obj.rhoList);

            obj.Cab = zeros(4,nA,nB,nRho);

            for iRho = 1:nRho

                rho = obj.rhoList(iRho);
                m0  = sqrt(1-rho);
                m1  = m0;
                m2  = m0;

                for iA = 1:nA

                    a = obj.paramA(iA);

                    for iB = 1:nB

                        b = obj.paramB(iB);

                        obj.createLattice(a,b);

                        Ch = obj.computeHomogenizationFromM(m1,m2);
                        obj.Cab(:,iA,iB,iRho) = obj.extractComponents(Ch);

                    end
                end
            end

        end

        function computeBetaCompetition(obj)

            betaMin = atan(obj.minParamB);
            betaMax = atan(obj.maxParamB);

            obj.paramBeta = sort(unique([ ...
                linspace(betaMin,betaMax,obj.nStepsBeta), ...
                pi/4]));

            nBeta = numel(obj.paramBeta);
            nRho  = numel(obj.rhoList);

            obj.Cbeta = zeros(4,nBeta,nRho);

            for iRho = 1:nRho

                rho = obj.rhoList(iRho);
                m0  = sqrt(1-rho);
                m1  = m0;
                m2  = m0;

                for iBeta = 1:nBeta

                    beta = obj.paramBeta(iBeta);
                    b    = tan(beta);

                    obj.createLattice(0,b);

                    Ch = obj.computeHomogenizationFromM(m1,m2);
                    obj.Cbeta(:,iBeta,iRho) = obj.extractComponents(Ch);

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

            names = {'C_{11}','C_{22}','C_{12}','C_{33}'};
            colors = lines(numel(obj.rhoList));
            leg = arrayfun(@(r) sprintf('\\rho = %.1f',r),obj.rhoList,'UniformOutput',false);

            figure( ...
                'Color','w', ...
                'Position',[100 100 1250 900]);

            t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

            for i = 1:4

                ax = nexttile;
                hold(ax,'on');

                h = gobjects(numel(obj.rhoList),1);

                for iRho = 1:numel(obj.rhoList)

                    xiVec = obj.paramXi{iRho};
                    data  = obj.Cxi{iRho}(i,:) / obj.Csolid(i);

                    h(iRho) = plot( ...
                        ax, ...
                        xiVec, ...
                        data, ...
                        'LineWidth',2, ...
                        'Color',colors(iRho,:));

                end

                xline(ax,1,':');
                yline(ax,1,':');

                xlabel(ax,'\xi');
                ylabel(ax,'C_{IJ}^{h}/C_{IJ}^{solid}');
                title(ax,names{i});
                grid(ax,'on');
                box(ax,'on');

                if i == 1
                    legend(ax,h,leg,'Location','best');
                end

            end

            title(t,'Dimensionless C_{IJ}(\xi), a = 0, b = 1');

        end

        function plotACompetition(obj)

            names = {'C_{11}','C_{22}','C_{12}','C_{33}'};
            colors = lines(numel(obj.rhoList));
            leg = arrayfun(@(r) sprintf('\\rho = %.1f',r),obj.rhoList,'UniformOutput',false);

            figure( ...
                'Color','w', ...
                'Position',[110 110 1250 900]);

            t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

            for i = 1:4

                ax = nexttile;
                hold(ax,'on');

                h = gobjects(numel(obj.rhoList),1);

                for iRho = 1:numel(obj.rhoList)

                    data = squeeze(obj.Ca(i,:,iRho)) / obj.Csolid(i);

                    h(iRho) = plot( ...
                        ax, ...
                        obj.paramA, ...
                        data, ...
                        'LineWidth',2, ...
                        'Color',colors(iRho,:));

                end

                xline(ax,0,':');
                yline(ax,1,':');

                xlabel(ax,'a');
                ylabel(ax,'C_{IJ}^{h}/C_{IJ}^{solid}');
                title(ax,names{i});
                grid(ax,'on');
                box(ax,'on');

                if i == 1
                    legend(ax,h,leg,'Location','best');
                end

            end

            title(t,'Dimensionless C_{IJ}(a), \xi = 1, b = 1');

        end

        function plotBCompetition(obj)

            names = {'C_{11}','C_{22}','C_{12}','C_{33}'};
            colors = lines(numel(obj.rhoList));
            leg = arrayfun(@(r) sprintf('\\rho = %.1f',r),obj.rhoList,'UniformOutput',false);

            figure( ...
                'Color','w', ...
                'Position',[120 120 1250 900]);

            t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

            for i = 1:4

                ax = nexttile;
                hold(ax,'on');

                h = gobjects(numel(obj.rhoList),1);

                for iRho = 1:numel(obj.rhoList)

                    data = squeeze(obj.Cb(i,:,iRho)) / obj.Csolid(i);

                    h(iRho) = plot( ...
                        ax, ...
                        obj.paramB, ...
                        data, ...
                        'LineWidth',2, ...
                        'Color',colors(iRho,:));

                end

                xline(ax,1,':');
                yline(ax,1,':');

                xlabel(ax,'b');
                ylabel(ax,'C_{IJ}^{h}/C_{IJ}^{solid}');
                title(ax,names{i});
                grid(ax,'on');
                box(ax,'on');

                if i == 1
                    legend(ax,h,leg,'Location','best');
                end

            end

            title(t,'Dimensionless C_{IJ}(b), a = 0, \xi = 1');

        end

        function plotABCompetition(obj)

            names = {'C_{11}','C_{22}','C_{12}','C_{33}'};
            colors = lines(numel(obj.rhoList));
            leg = arrayfun(@(r) sprintf('\\rho = %.1f',r),obj.rhoList,'UniformOutput',false);

            [B,A] = meshgrid(obj.paramB,obj.paramA);

            figure( ...
                'Color','w', ...
                'Position',[130 80 1300 950]);

            t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

            for i = 1:4

                ax = nexttile;
                hold(ax,'on');

                h = gobjects(numel(obj.rhoList),1);

                for iRho = 1:numel(obj.rhoList)

                    data = squeeze(obj.Cab(i,:,:,iRho)) / obj.Csolid(i);

                    h(iRho) = mesh( ...
                        ax, ...
                        A, ...
                        B, ...
                        data, ...
                        'EdgeColor',colors(iRho,:), ...
                        'FaceColor','none', ...
                        'LineWidth',1.2);

                end

                xlabel(ax,'a');
                ylabel(ax,'b');
                zlabel(ax,'C_{IJ}^{h}/C_{IJ}^{solid}');
                title(ax,names{i});
                grid(ax,'on');
                box(ax,'on');
                view(ax,45,30);

                if i == 1
                    legend(ax,h,leg,'Location','best');
                end

            end

            title(t,'Dimensionless C_{IJ}(a,b), \xi = 1');

        end

        function plotBetaCompetition(obj)

            names = {'C_{11}','C_{22}','C_{12}','C_{33}'};
            colors = lines(numel(obj.rhoList));
            leg = arrayfun(@(r) sprintf('\\rho = %.1f',r),obj.rhoList,'UniformOutput',false);

            figure( ...
                'Color','w', ...
                'Position',[140 140 1250 900]);

            t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

            for i = 1:4

                ax = nexttile;
                hold(ax,'on');

                h = gobjects(numel(obj.rhoList),1);

                for iRho = 1:numel(obj.rhoList)

                    data = squeeze(obj.Cbeta(i,:,iRho)) / obj.Csolid(i);

                    h(iRho) = plot( ...
                        ax, ...
                        obj.paramBeta, ...
                        data, ...
                        'LineWidth',2, ...
                        'Color',colors(iRho,:));

                end

                xline(ax,pi/4,':');
                yline(ax,1,':');

                xlabel(ax,'\beta');
                ylabel(ax,'C_{IJ}^{h}/C_{IJ}^{solid}');
                title(ax,names{i});
                grid(ax,'on');
                box(ax,'on');

                if i == 1
                    legend(ax,h,leg,'Location','best');
                end

            end

            title(t,'Dimensionless C_{IJ}(\beta), a = 0, \xi = 1, b = tan(\beta)');

        end

    end

end