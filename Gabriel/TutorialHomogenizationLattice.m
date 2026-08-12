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
        maxParamB
        maxParamRho
    end

    methods (Access = public)

        function obj = TutorialHomogenizationLattice()
            obj.init();
            obj.computeHoleParams();
            obj.compute();
            obj.fitting();
            obj.plot();
            
        end

    end

    methods (Access = private)

        function init(obj)
            obj.E           = 1;
            obj.nu          = 0.3;
            obj.meshType    = 'Square';
            obj.meshN       = 80;
            obj.holeType    = 'Square';
            obj.pnorm       = 'Inf';
            obj.nStepsB     = 60;
            obj.nStepsRho   = 60;
            obj.monitoring  = false;
            obj.maxParamB   = 0.6;
            obj.maxParamRho = 0.979;
        end

        function computeHoleParams(obj)
            obj.paramB   = linspace(-0.6,    obj.maxParamB,   obj.nStepsB);
            obj.paramRho = linspace(1e-9, obj.maxParamRho, obj.nStepsRho);
        end

        function compute(obj)
            nB   = length(obj.paramB);
            nRho = length(obj.paramRho);
            mat  = zeros(2, 2, 2, 2, nRho, nB);
            volF = zeros(nRho, nB);

            for iRho = 1:nRho
                rho_val = obj.paramRho(iRho);
                obj.currentRho = rho_val;
                fprintf('\n=== rho = %.4f ===\n', rho_val);

                for iB = 1:nB
                    b_val = obj.paramB(iB);
                    obj.currentB = b_val;

                   
                    a  = exp(b_val^2);
                    d  = (1 + b_val^2) / a;
                    v1 = [a,     b_val];
                    v2 = [b_val, d    ];
                    obj.latticeVectors = [v1; v2];

                    obj.defineMesh();

                    mat(:,:,:,:,iRho,iB) = obj.computeHomogenization(rho_val, b_val);
                    volF(iRho,iB)        = obj.computeVolumeFraction(rho_val, b_val);

                    if mod(iB, 5) == 0 || iB == nB
                        fprintf('  b = %.4f  volF = %.4f\n', b_val, volF(iRho,iB));
                    end
                end
            end

            obj.Chomog  = mat;
            obj.volFrac = volF;
        end

        function matHomog = computeHomogenization(obj, rho_val, b_val)
            dens     = obj.createDensityLevelSet(rho_val, b_val);
            mat      = obj.createDensityMaterial(dens);
            matHomog = obj.solveElasticMicroProblem(mat, dens);
        end

        function lsf = createDensityLevelSet(obj, rho_val,b_val)
           
            ls = obj.computeLevelSet(obj.baseMesh, rho_val);

            sUm.backgroundMesh = obj.baseMesh;
            sUm.boundaryMesh   = obj.baseMesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);

            ls  = CharacteristicFunction.create(uMesh);
            s.trial = obj.test;
            s.mesh  = obj.baseMesh;
            f   = FilterLump(s);
            lsf = f.compute(ls, 2);
        end

        function ls = computeLevelSet(obj, mesh, rho_val)
            
            gPar.type  = obj.holeType;
            gPar.pnorm = obj.pnorm;

            coord    = mesh.coord;
            center_x = (min(coord(:,1)) + max(coord(:,1))) / 2;
            center_y = (min(coord(:,2)) + max(coord(:,2))) / 2;
            gPar.xCoorCenter = center_x;
            gPar.yCoorCenter = center_y;

            v1  = obj.latticeVectors(1, :);
            v2  = obj.latticeVectors(2, :);
            phi = atan2(v1(2), v1(1));
            gPar.a1       = v1;
            gPar.a2       = v2;
            gPar.rotation = phi;

            switch obj.holeType
                case 'Circle'
                    gPar.radius = rho_val / 2;
                case 'Square'
                    gPar.length = sqrt(1 - rho_val);
                case 'SmoothRectangle'
                    gPar.xSide = rho_val;
                    gPar.ySide = rho_val / 2;
                    gPar.pnorm = 16;
                case 'Ellipse'
                    gPar.type  = 'SmoothRectangle';
                    gPar.xSide = rho_val(1);
                    gPar.ySide = rho_val(2);
                    gPar.pnorm = 2;
                case 'SmoothHexagon'
                    gPar.radius = rho_val;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
                case 'ReinforcedHoneycomb'
                    gPar.theta    = 1 - rho_val;
                    gPar.eps      = 1;
                    gPar.normal   = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
                    gPar.radius   = rho_val;
                    gPar.rotation = phi;
            end

            g      = GeometricalFunction(gPar);
            phiFun = g.computeLevelSetFunction(mesh);
            ls     = -phiFun.fValues;
        end

        function defineMesh(obj)
            s.latticeVectors = obj.latticeVectors;
            s.divUnit        = obj.meshN;
            s.filename       = '';
            MC = MeshCreator(s);
            MC.computeMeshNodes();
            s.coord         = MC.coord;
            s.connec        = MC.connec;
            obj.baseMesh    = Mesh.create(s);
            obj.masterSlave = MC.masterSlaveIndex;
            
            obj.test        = LagrangianFunction.create(obj.baseMesh, 1, 'P1');
            obj.Mmass       = IntegrateLHS(@(u,v) DP(v,u), obj.test, obj.test, obj.baseMesh, 'Domain');
        end

        function mat = createDensityMaterial(obj, lsf)
            s.interpolation  = 'SIMPALL';
            s.dim            = '2D';
            s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1e-6*obj.E, obj.nu, obj.baseMesh.ndim);
            s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1e-6*obj.E, obj.nu);
            s.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(obj.E, obj.nu, obj.baseMesh.ndim);
            s.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(obj.E, obj.nu);
            mI = MaterialInterpolator.create(s);

            x{1} = lsf;
            s.mesh                 = obj.baseMesh;
            s.type                 = 'DensityBased';
            s.density              = x;
            s.materialInterpolator = mI;
            s.dim                  = '2D';
            mat = Material.create(s);
        end

        function matHomog = solveElasticMicroProblem(obj, material, dens)
            if obj.monitoring
                close all
                dens.plot
                shading interp
                colormap(flipud(pink))
                drawnow
            end

            s.mesh               = obj.baseMesh;
            s.material           = material;
            s.scale              = 'MICRO';
            s.dim                = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(obj.baseMesh);
            s.solverCase         = DirectSolver();
            s.solverType         = 'REDUCED';
            s.solverMode         = 'FLUC';
            fem = ElasticProblemMicro(s);
            material.setDesignVariable({dens})
            fem.updateMaterial(material.obtainTensor())
            fem.solve();

            totVol   = obj.baseMesh.computeVolume();
            matHomog = fem.Chomog / totVol;
        end

        function bc = createBoundaryConditions(obj, mesh)

            coord = mesh.coord;

            cornerCoord = coord(1:4,:);

            tol2 = 1e-20;

            isCorner = @(coor) ...
                (sum((coor - cornerCoord(1,:)).^2,2) < tol2) | ...
                (sum((coor - cornerCoord(2,:)).^2,2) < tol2) | ...
                (sum((coor - cornerCoord(3,:)).^2,2) < tol2) | ...
                (sum((coor - cornerCoord(4,:)).^2,2) < tol2);

            sDir{1}.domain    = @(coor) isCorner(coor);
            sDir{1}.direction = [1, 2];
            sDir{1}.value     = 0;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dirichletFun = [dirichletFun, ...
                    DirichletCondition(mesh, sDir{i})];
            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = 1;
            s.mesh         = mesh;

            bc = BoundaryConditions(s);
            bc.updatePeriodicConditions(obj.masterSlave);

        end

        function fracVol = computeVolumeFraction(obj, rho_val, b_val)
            rho    = obj.createDensityLevelSet(rho_val, b_val);
            volDom = Integrator.compute(ConstantFunction.create(1, obj.baseMesh), obj.baseMesh, 2);
            fracVol = Integrator.compute(rho, rho.mesh, 2) / volDom;
        end

        % function fitting(obj)
        %     [obj.f, obj.df, obj.ddf] = DamageHomogenizationFitter.computeNN( ...
        %         obj.paramB, obj.paramRho, obj.Chomog);
        % end
        function fitting(obj)

            s.retrain      = true;   
            s.pol_deg      = 6;      
            s.hiddenLayers = [150 200 300 200 150 50];
            % s.hiddenLayers = [50 100 50];
            s.maxEpochs    = 300000;
            s.learningRate = 0.02;

            [obj.f, obj.df, ~] = DamageHomogenizationFitter.computeNN(obj.paramB, obj.paramRho, obj.Chomog, s);

        end

        function plot(obj)
            [B, R] = meshgrid(obj.paramB, obj.paramRho);
            components = {[1,1,1,1], 'C_{1111}'; ...
                          [2,2,2,2], 'C_{2222}'; ...
                          [1,2,1,2], 'C_{1212}'};
            figure;
            tiledlayout(1, 3, 'TileSpacing', 'compact');
            for k = 1:3
                idx  = components{k,1};
                data = squeeze(obj.Chomog(idx(1),idx(2),idx(3),idx(4),:,:));
                nexttile
                surf(B, R, data, 'EdgeColor', 'none')
                xlabel('b'); ylabel('\rho'); zlabel(components{k,2})
                title(components{k,2}); grid on; view(45, 30)
            end
        end

        

    end

end