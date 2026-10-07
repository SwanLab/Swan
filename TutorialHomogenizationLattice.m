classdef TutorialHomogenizationLattice < handle

    properties (Access = public)
        paramB
        paramRho
        Chomog
        volFrac
        f, df, ddf
    end

    properties (Access = private)
        E, nu
        meshType
        meshN
        holeType
        nStepsB
        nStepsRho
        pnorm
        monitoring
        Mmass
        currentB
        currentRho
        latticeVectors
        baseMesh
        masterSlave
        test
        minParamB
        maxParamB
        maxParamRho
        anisotropyCoupling
    end

    methods (Access = public)

        function obj = TutorialHomogenizationLattice()
            obj.init();
            obj.computeHoleParams();
            obj.compute();
            obj.fitting();
            obj.plot();
        end

        function exportMicrosCurrentB(obj,rhoVal,b_vals)

            if nargin < 3
                b_vals = [0.3 0.6 1 1.5 2 3];
            end

            if nargin < 2
                rhoVal = 0.5;
            end

            folderName = sprintf('Micros_B_rho_%0.2f',rhoVal);

            if ~exist(folderName,'dir')
                mkdir(folderName);
            end

            a_fixed = 0;
            margin = 0.03;

            for i = 1:numel(b_vals)

                bVal = b_vals(i);

                k = 1/sqrt(1-a_fixed^2);

                B = k * [1,       a_fixed;
                         a_fixed, 1];

                D = [bVal,     0;
                     0,   1/bVal];

                F = B * D;

                v1 = F(1,:);
                v2 = F(2,:);

                obj.latticeVectors = [v1;v2];

                obj.defineMesh();

                dens = obj.createDensityLevelSet(rhoVal,bVal);

                funP0   = dens.project('P0');
                rhoElem = squeeze(funP0.fValues);

                coordPlot = obj.baseMesh.coord;

                xc = 0.5*(min(coordPlot(:,1)) + max(coordPlot(:,1)));
                yc = 0.5*(min(coordPlot(:,2)) + max(coordPlot(:,2)));

                coordPlot(:,1) = coordPlot(:,1) - xc;
                coordPlot(:,2) = coordPlot(:,2) - yc;

                xMin = min(coordPlot(:,1));
                xMax = max(coordPlot(:,1));
                yMin = min(coordPlot(:,2));
                yMax = max(coordPlot(:,2));

                dx = xMax - xMin;
                dy = yMax - yMin;

                halfRange = 0.5*max(dx,dy)*(1 + 2*margin);

                fig = figure( ...
                    'Color','w', ...
                    'Position',[100 100 900 900]);

                ax = axes(fig);
                hold(ax,'on');

                patch(ax, ...
                    'Faces',obj.baseMesh.connec, ...
                    'Vertices',coordPlot, ...
                    'FaceVertexCData',rhoElem, ...
                    'FaceColor','flat', ...
                    'EdgeColor','none');

                colormap(ax,flipud(gray(256)));
                caxis(ax,[0 1]);

                axis(ax,'equal');

                xlim(ax,[-halfRange halfRange]);
                ylim(ax,[-halfRange halfRange]);

                grid(ax,'on');
                box(ax,'on');

                xlabel(ax,'x');
                ylabel(ax,'y');

                title(ax, ...
                    sprintf('$b = %.2f,\\; \\rho = %.2f$', ...
                    bVal,rhoVal), ...
                    'Interpreter','latex', ...
                    'FontSize',18);

                set(ax,'FontSize',13,'Layer','top');

                bString = sprintf('%0.2f',bVal);
                bString = strrep(bString,'.','p');

                rhoString = sprintf('%0.2f',rhoVal);
                rhoString = strrep(rhoString,'.','p');

                fileName = fullfile(folderName, ...
                    sprintf('Micro_b_%s_rho_%s.png', ...
                    bString,rhoString));

                exportgraphics( ...
                    fig,fileName, ...
                    'Resolution',300);

                close(fig);

            end
        end

    end

    methods (Access = private)

        function init(obj)

            obj.E           = 1;
            obj.nu          = 0.3;
            obj.meshType    = 'Square';
            obj.meshN       = 60;
            obj.holeType    = 'RectangleAffine';
            obj.pnorm       = 'Inf';

            obj.nStepsB     = 73;
            obj.nStepsRho   = 71;

            obj.monitoring  = false;

            obj.maxParamB = 3.0;
            obj.minParamB = 1/obj.maxParamB;

            obj.maxParamRho = 0.998;

            obj.anisotropyCoupling = 0.9;

        end

        function computeHoleParams(obj)

            obj.paramB = linspace( ...
                obj.minParamB, ...
                obj.maxParamB, ...
                obj.nStepsB);

            obj.paramRho = linspace( ...
                1e-9, ...
                obj.maxParamRho, ...
                obj.nStepsRho);

        end

        function compute(obj)

            nB   = length(obj.paramB);
            nRho = length(obj.paramRho);

            mat  = zeros(2,2,2,2,nRho,nB);
            volF = zeros(nRho,nB);

            a_fixed = 0;

            for iRho = 1:nRho

                rho_val = obj.paramRho(iRho);
                obj.currentRho = rho_val;

                fprintf('\n=== rho = %.4f ===\n',rho_val);

                for iB = 1:nB

                    b_val = obj.paramB(iB);
                    obj.currentB = b_val;

                    if b_val <= 0
                        error('Parameter b must satisfy b > 0.');
                    end

                    k = 1/sqrt(1-a_fixed^2);

                    B = k * [1,       a_fixed;
                             a_fixed, 1];

                    D = [b_val,       0;
                         0,      1/b_val];

                    F = B * D;

                    v1 = F(1,:);
                    v2 = F(2,:);

                    obj.latticeVectors = [v1;
                                          v2];

                    obj.defineMesh();

                    mat(:,:,:,:,iRho,iB) = ...
                        obj.computeHomogenization( ...
                        rho_val,b_val);

                    volF(iRho,iB) = ...
                        obj.computeVolumeFraction( ...
                        rho_val,b_val);

                    if mod(iB,5) == 0 || iB == nB

                        fprintf([ ...
                            '  a = %.1f   b = %.4f   ' ...
                            'det(F) = %.8f   volF = %.4f\n'], ...
                            a_fixed, ...
                            b_val, ...
                            det(F), ...
                            volF(iRho,iB));

                    end

                end
            end

            obj.Chomog  = mat;
            obj.volFrac = volF;

        end

        function matHomog = ...
                computeHomogenization(obj,rho_val,b_val)

            dens = ...
                obj.createDensityLevelSet( ...
                rho_val,b_val);

            mat = ...
                obj.createDensityMaterial(dens);

            matHomog = ...
                obj.solveElasticMicroProblem( ...
                mat,dens);

        end

        function lsf = ...
                createDensityLevelSet(obj,rho_val,b_val)

            ls = ...
                obj.computeLevelSet( ...
                obj.baseMesh, ...
                rho_val, ...
                b_val);

            sUm.backgroundMesh = obj.baseMesh;

            sUm.boundaryMesh = ...
                obj.baseMesh.createBoundaryMesh;

            uMesh = UnfittedMesh(sUm);

            uMesh.compute(ls);

            ls = CharacteristicFunction.create(uMesh);

            s.trial = obj.test;
            s.mesh  = obj.baseMesh;

            f = FilterLump(s);

            lsf = f.compute(ls,2);

        end

        function ls = ...
                computeLevelSet(obj,mesh,rho_val,b_val)

            if rho_val < 0 || rho_val >= 1

                error( ...
                    'rho must satisfy 0 <= rho < 1.');

            end

            if b_val <= 0

                error( ...
                    'b must satisfy b > 0.');

            end

            gPar.type = obj.holeType;

            coord = mesh.coord;

            center_x = ...
                (min(coord(:,1)) + ...
                max(coord(:,1)))/2;

            center_y = ...
                (min(coord(:,2)) + ...
                max(coord(:,2)))/2;

            gPar.xCoorCenter = center_x;
            gPar.yCoorCenter = center_y;

            v1 = obj.latticeVectors(1,:);
            v2 = obj.latticeVectors(2,:);

            gPar.a1 = v1;
            gPar.a2 = v2;

            eta = ...
                obj.anisotropyCoupling * ...
                (b_val - 1/b_val) / ...
                (obj.maxParamB - 1/obj.maxParamB);

            eta = max( ...
                -obj.anisotropyCoupling, ...
                min(obj.anisotropyCoupling,eta));

            voidFraction = 1-rho_val;

            m1 = ...
                voidFraction.^((1-eta)/2);

            m2 = ...
                voidFraction.^((1+eta)/2);

            switch obj.holeType

                case 'RectangleAffine'

                    gPar.xSide = m1;
                    gPar.ySide = m2;

                case 'Square'

                    gPar.length = ...
                        sqrt(1-rho_val);

                case 'Circle'

                    gPar.radius = rho_val/2;

                case 'CrossedBars'

                    gPar.width = ...
                        0.25*( ...
                        1-sqrt(max(0,1-rho_val)));

                case 'TwoHorizontalBars'

                    rhoSat = 0.15;
                    tFrameMax = 0.02;

                    if rho_val <= rhoSat

                        tFrame = ...
                            tFrameMax * ...
                            rho_val/rhoSat;

                    else

                        tFrame = tFrameMax;

                    end

                    AFrame = ...
                        1-(1-2*tFrame)^2;

                    remaining = ...
                        max(rho_val-AFrame,0);

                    tBar = ...
                        remaining / ...
                        (2*(1-2*tFrame));

                    tBar = ...
                        min( ...
                        tBar, ...
                        0.5-tFrame);

                    gPar.width = tBar;

                    gPar.frameWidth = ...
                        tFrame;

                case 'SmoothRectangle'

                    gPar.xSide = rho_val;
                    gPar.ySide = rho_val/2;
                    gPar.pnorm = 16;

                    phi = atan2( ...
                        v1(2), ...
                        v1(1));

                    gPar.rotation = phi;

                case 'Ellipse'

                    gPar.type = ...
                        'SmoothRectangle';

                    gPar.xSide = rho_val(1);
                    gPar.ySide = rho_val(2);

                    gPar.pnorm = 2;

                    phi = atan2( ...
                        v1(2), ...
                        v1(1));

                    gPar.rotation = phi;

                case 'SmoothHexagon'

                    gPar.radius = rho_val;

                    gPar.normal = [ ...
                        0 1;
                        sqrt(3)/2 1/2;
                        sqrt(3)/2 -1/2];

                case 'ReinforcedHoneycomb'

                    gPar.theta = ...
                        1-rho_val;

                    gPar.eps = 1;

                    gPar.normal = [ ...
                        0 1;
                        sqrt(3)/2 1/2;
                        sqrt(3)/2 -1/2];

                    gPar.radius = rho_val;

                    phi = atan2( ...
                        v1(2), ...
                        v1(1));

                    gPar.rotation = phi;

            end

            g = GeometricalFunction(gPar);

            phiFun = ...
                g.computeLevelSetFunction(mesh);

            switch obj.holeType

                case { ...
                        'CrossedBars', ...
                        'TwoHorizontalBars'}

                    ls = phiFun.fValues;

                otherwise

                    ls = -phiFun.fValues;

            end

        end

        function defineMesh(obj)

            s.latticeVectors = ...
                obj.latticeVectors;

            s.divUnit  = obj.meshN;
            s.filename = '';

            MC = MeshCreator(s);

            MC.computeMeshNodes();

            s.coord  = MC.coord;
            s.connec = MC.connec;

            obj.baseMesh = Mesh.create(s);

            obj.masterSlave = ...
                MC.masterSlaveIndex;

            obj.test = ...
                LagrangianFunction.create( ...
                obj.baseMesh,1,'P1');

            obj.Mmass = ...
                IntegrateLHS( ...
                @(u,v) DP(v,u), ...
                obj.test, ...
                obj.test, ...
                obj.baseMesh, ...
                'Domain');

        end

        function mat = createDensityMaterial(obj,lsf)

            s.interpolation = 'SIMPALL';
            s.dim           = '2D';

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

        function matHomog = ...
                solveElasticMicroProblem(obj,material,dens)

            if obj.monitoring

                close all

                dens.plot

                shading interp

                colormap(flipud(pink))

                drawnow

            end

            s.mesh     = obj.baseMesh;
            s.material = material;
            s.scale    = 'MICRO';
            s.dim      = '2D';

            s.boundaryConditions = ...
                obj.createBoundaryConditions( ...
                obj.baseMesh);

            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';

            fem = ElasticProblemMicro(s);

            material.setDesignVariable({dens});

            fem.updateMaterial( ...
                material.obtainTensor());

            fem.solve();

            totVol = ...
                obj.baseMesh.computeVolume();

            matHomog = ...
                fem.Chomog/totVol;

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

            sDir{1}.domain = ...
                @(coor) isCorner(coor);

            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            dirichletFun = [];

            for i = 1:numel(sDir)

                dirichletFun = ...
                    [dirichletFun, ...
                    DirichletCondition( ...
                    mesh,sDir{i})];

            end

            s.dirichletFun = dirichletFun;

            s.pointloadFun = [];
            s.periodicFun  = 1;
            s.mesh         = mesh;

            bc = BoundaryConditions(s);

            bc.updatePeriodicConditions( ...
                obj.masterSlave);

        end

        function fracVol = ...
                computeVolumeFraction(obj,rho_val,b_val)

            rho = ...
                obj.createDensityLevelSet( ...
                rho_val,b_val);

            volDom = ...
                Integrator.compute( ...
                ConstantFunction.create( ...
                1,obj.baseMesh), ...
                obj.baseMesh,2);

            fracVol = ...
                Integrator.compute( ...
                rho,rho.mesh,2)/volDom;

        end

        function fitting(obj)

            Cfit = permute( ...
                obj.Chomog, ...
                [1 2 3 4 6 5]);

            paramVectors = { ...
                obj.paramB, ...
                obj.paramRho};

            s.retrain = true;

            s.parameterNames = { ...
                'b', ...
                'rho'};

            s.transforms = { ...
                'identity', ...
                'identity'};

            s.featureMap = 'direct';

            s.hiddenLayers = ...
                [150 200 300 200 150 50];

            s.HUtype = 'tanh';
            s.OUtype = 'linear';

            s.maxEpochs = 100000;
            s.learningRate = 0.01;

            s.lambda = 0;
            s.seed = 0;

            s.saveFile = ...
                'HomogNN_CaseC_brho_identity_tanh_100k.mat';

            s.historyFile = ...
                'NNhistory_CaseC_brho_identity_tanh_100k.mat';

            s.referenceDataFile = '';

            s.allowExtrapolation = false;

            [obj.f,obj.df,~] = ...
                DamageHomogenizationFitter.computeNN( ...
                paramVectors, ...
                Cfit, ...
                s);

        end

        function plot(obj)

            [B,R] = ...
                meshgrid( ...
                obj.paramB, ...
                obj.paramRho);

            components = { ...
                [1,1,1,1], 'C_{1111}'; ...
                [2,2,2,2], 'C_{2222}'; ...
                [1,2,1,2], 'C_{1212}'};

            figure;

            tiledlayout( ...
                1,3, ...
                'TileSpacing','compact');

            for k = 1:3

                idx = components{k,1};

                data = squeeze( ...
                    obj.Chomog( ...
                    idx(1),idx(2), ...
                    idx(3),idx(4),:,:));

                nexttile

                surf( ...
                    B,R,data, ...
                    'EdgeColor','none');

                xlabel('b');
                ylabel('\rho');
                zlabel(components{k,2});

                title(components{k,2});

                grid on

                view(45,30);

            end

        end

    end

end