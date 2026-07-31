classdef TutorialHomogenizationLatticeBOnly < handle
    % Caso 1 variavel: rho fixo em 0.5, b varia.
    % Gera vademecum 1D (so em b) e ajusta com polinomio (computePolynomial).

    properties (Access = public)
        paramB
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
        pnorm
        monitoring
        Mmass
        currentB
        fixedRho
        latticeVectors
        baseMesh
        masterSlave
        test
        maxParamB
        degPoly
    end

    methods (Access = public)

        function obj = TutorialHomogenizationLatticeBOnly()
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
            obj.meshN       = 100;
            obj.holeType    = 'Square';
            obj.pnorm       = 'Inf';
            obj.nStepsB     = 81;
            obj.monitoring  = false;
            obj.maxParamB   = 0.6;
            obj.fixedRho    = 0.4;    
            obj.degPoly     = 6;
        end

        function computeHoleParams(obj)
            obj.paramB = linspace(0, obj.maxParamB, obj.nStepsB);
        end

        function compute(obj)

            nB   = length(obj.paramB);
            mat  = zeros(2,2,2,2,nB);
            volF = zeros(nB,1);

            rho_val = obj.fixedRho;

            

            for iB = 1:nB

                
                beta_val = obj.paramB(iB);

                
                b_geom = obj.maxParamB-abs(beta_val);

                obj.currentB = b_geom;

                a = exp(b_geom^2);
                d = (1+b_geom^2)/a;

                v1 = [a,      b_geom];
                v2 = [b_geom, d     ];

                obj.latticeVectors = [v1;v2];

                obj.defineMesh();

                mat(:,:,:,:,iB) = ...
                    obj.computeHomogenization(rho_val,b_geom);

                volF(iB) = ...
                    obj.computeVolumeFraction(rho_val,b_geom);

                if mod(iB,5) == 0 || iB == 1 || iB == nB

                    fprintf([ ...
                        'beta = %.4f  ->  bGeom = %.4f', ...
                        '  volF = %.4f\n'], ...
                        beta_val,b_geom,volF(iB));

                end

            end

            obj.Chomog  = mat;
            obj.volFrac = volF;

            fprintf('============================================================\n');

        end

        function matHomog = computeHomogenization(obj, rho_val, b_val)
            dens     = obj.createDensityLevelSet(rho_val, b_val);
            mat      = obj.createDensityMaterial(dens);
            matHomog = obj.solveElasticMicroProblem(mat, dens);
        end

        function lsf = createDensityLevelSet(obj, rho_val, b_val)
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
                case 'Square'
                    gPar.length = sqrt(1 - rho_val);
                case 'Circle'
                    gPar.radius = rho_val / 2;
                otherwise
                    error('holeType nao suportado neste tutorial: %s', obj.holeType);
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

            obj.test  = LagrangianFunction.create(obj.baseMesh, 1, 'P1');
            obj.Mmass = IntegrateLHS(@(u,v) DP(v,u), obj.test, obj.test, obj.baseMesh, 'Domain');
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
            v1 = obj.latticeVectors(1,:);
            v2 = obj.latticeVectors(2,:);
            origin = [min(coord(:,1)), min(coord(:,2))];

            [~,iFix] = min(sum((coord - origin).^2,2));
            isFix = @(coor) sum((coor - coord(iFix,:)).^2,2) < 1e-20;

            sDir{1}.domain    = @(coor) isFix(coor);
            sDir{1}.direction = [1, 2];
            sDir{1}.value     = 0;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dirichletFun = [dirichletFun, DirichletCondition(mesh, sDir{i})];
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

        function fitting(obj)
            % Fit polinomial 1D (so em b), reaproveitando o metodo
            % computePolynomial ja existente no DamageHomogenizationFitter
            [obj.f, obj.df, obj.ddf] = DamageHomogenizationFitter.computePolynomial( ...
                obj.degPoly, obj.paramB, obj.Chomog);
        end

        function plot(obj)
            components = {[1,1,1,1], 'C_{1111}'; ...
                          [2,2,2,2], 'C_{2222}'; ...
                          [1,2,1,2], 'C_{1212}'};
            figure;
            tiledlayout(1, 3, 'TileSpacing', 'compact');
            for k = 1:3
                idx  = components{k,1};
                data = squeeze(obj.Chomog(idx(1),idx(2),idx(3),idx(4),:));
                nexttile
                plot(obj.paramB, data, 'LineWidth', 1.5)
                xlabel('b'); ylabel(components{k,2})
                title(components{k,2}); grid on
            end
        end

    end

end