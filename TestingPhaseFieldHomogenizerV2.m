classdef TestingPhaseFieldHomogenizerV2 < handle

    properties (Access = private)
        material
        mesh
        hole
    end

    properties (Access = private)
        baseMesh
        masterSlave
        bc
        matInterp
        densTestFun
        filter
        monitor
    end

    methods (Access = public)
        
        function obj = TestingPhaseFieldHomogenizerV2(cParams)
            obj.init(cParams)
            obj.defineMesh();
            obj.defineBoundaryConditions();
            obj.defineMaterialInterpolation();
            obj.defineDensityFilter();
            obj.createMonitoring(cParams);
        end
        
        function [mat,comb] = compute(obj)
            [comb,nComb] = obj.computeHoleSizeCombinations();
            matComb = zeros(2,2,2,2,nComb);
            for i=1:nComb
                holeSize = comb(i,:);
                matComb(:,:,:,:,i) = obj.computeHomogenization(holeSize,i);
            end
            mat = reshape(matComb, [2,2,2,2,cellfun(@numel,obj.hole.size)]);  
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.material = cParams.material;
            obj.hole     = cParams.hole;
        end

        function defineMesh(obj)
            switch obj.mesh.type
                case 'Square'
                    obj.baseMesh.center = [0.5,0.5];
                    s.c = [1,1];
                    s.theta = [0,90];
                    s.divUnit = obj.mesh.N;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
                case 'Hexagon'
                    obj.baseMesh.center = [0.5,sqrt(1-0.5^2)];
                    s.c = [1,1,1];
                    s.theta = [0,60,120];
                    s.divUnit = obj.mesh.N;
                    s.filename = '';
                    MC = MeshCreator(s);
                    MC.computeMeshNodes();
            end
            s.coord  = MC.coord;
            s.connec = MC.connec;
            obj.baseMesh.mesh = Mesh.create(s);
            obj.masterSlave   = MC.masterSlaveIndex;
        end

        function defineBoundaryConditions(obj)
            switch obj.mesh.type
                case 'Square'
                    isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
                    isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12);
                    isVertex = @(coor) (isTop(coor) & isLeft(coor));
                    sDir.domain    = @(coor) isVertex(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
                case 'Hexagon'
                    isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12);
                    sDir.domain    = @(coor) isLeft(coor);
                    sDir.direction = [1,2];
                    sDir.value     = 0;
            end
            s.mesh = obj.baseMesh.mesh;
            s.dirichletFun = DirichletCondition(obj.baseMesh.mesh,sDir);
            s.pointloadFun = [];
            s.periodicFun  = [];
            obj.bc = BoundaryConditions(s);
            obj.bc.updatePeriodicConditions(obj.masterSlave);
        end

        function defineMaterialInterpolation(obj)
            ndim = obj.baseMesh.mesh.ndim;
            matA = obj.material.A; 
            matB = obj.material.B;
            s.interpolation  = 'SIMPALL';
            s.dim            = ndim;
            s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(matA.E,matA.nu,ndim);
            s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(matA.E,matA.nu);
            s.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(matB.E,matB.nu,ndim);
            s.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(matB.E,matB.nu);
            obj.matInterp = MaterialInterpolator.create(s);
        end

        function defineDensityFilter(obj)
            obj.densTestFun = LagrangianFunction.create(obj.baseMesh.mesh,1,'P1');
            s.trial = obj.densTestFun;
            s.mesh  = obj.baseMesh.mesh;
            obj.filter = FilterLump(s);
        end

        function createMonitoring(obj,params)
            if params.monitoring
                s.shallDisplay = params.monitoring;
                s.funs = {obj.densTestFun};
                s.barLims = {[0,1]};
                s.maxNColumns = 1;
                s.titles = {'Density'};
                s.chartTypes = {'surf'};
                obj.monitor = Monitoring(s);
            else
                obj.monitor = [];
            end
        end

        function [comb,nComb] = computeHoleSizeCombinations(obj)
            nParams = numel(obj.hole.size);
            grids = cell(1,nParams);
            [grids{:}] = ndgrid(obj.hole.size{:});
            comb = zeros(numel(grids{1}),nParams);
            for i = 1:nParams
                comb(:,i) = grids{i}(:);
            end
            nComb = size(comb,1);
        end

        function matHomog = computeHomogenization(obj,l,step)
            density  = obj.createDensityLevelSet(l);
            mat      = obj.createDensityMaterial(density);
            matHomog = obj.solveElasticMicroProblem(mat,density);
            obj.updateMonitoring(density,step);
        end

        function lsf = createDensityLevelSet(obj,l)
            ls = obj.computeLevelSet(l);
            sUm.backgroundMesh = obj.baseMesh.mesh;
            sUm.boundaryMesh   = obj.baseMesh.mesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);

            ls = CharacteristicFunction.create(uMesh);
            lsf = obj.filter.compute(ls,2);
        end

        function ls = computeLevelSet(obj,size)
            gPar.type = obj.hole.type;
            gPar.xCoorCenter = obj.baseMesh.center(1);
            gPar.yCoorCenter = obj.baseMesh.center(2);
            switch obj.hole.type
                case 'Circle'
                    gPar.radius = size/2;
                case 'Square'
                    gPar.length = size;
                case 'Rectangle'
                    gPar.xSide  = size(1);
                    gPar.ySide  = size(2);
                case 'Crack'
                    gPar.type = 'Rectangle';
                    gPar.xSide = size;
                    gPar.ySide = (4/obj.mesh.N);
                case 'CircleCrack'
                    transSize = obj.hole.transSize;
                    sizeY = (size<transSize)*size + (size>=transSize)*transSize;
                    gPar.type = "SmoothRectangle";
                    gPar.xSide  = size;
                    gPar.ySide  = sizeY;
                    gPar.pnorm  = 2;
                case 'Ellipse'
                    gPar.type = "SmoothRectangle";
                    gPar.xSide  = size(1);
                    gPar.ySide  = size(2);
                    gPar.pnorm  = 2;
                case {'SmoothHexagon','Hexagon','ReinforcedHoneycomb'}
                    gPar.radius = size;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
                    gPar.pnorm = obj.hole.pnorm;
            end
            g   = GeometricalFunction(gPar);
            lsf = g.computeLevelSetFunction(obj.baseMesh.mesh);
            ls  = -lsf.fValues;
        end

        function mat = createDensityMaterial(obj,lsf)
            s.mesh                 = obj.baseMesh.mesh;
            s.type                 = 'DensityBased';
            s.density              = {lsf};
            s.materialInterpolator = obj.matInterp;
            s.dim                  = '2D';
            mat = Material.create(s);
        end

        function matHomog = solveElasticMicroProblem(obj,mat,dens)
            s.mesh = obj.baseMesh.mesh;
            s.material = mat;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.bc;
            s.solverCase = DirectSolver();
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);

            mat.setDesignVariable({dens})
            fem.updateMaterial(mat.obtainTensor())
            fem.solve();

            totVol = obj.baseMesh.mesh.computeVolume();
            matHomog = fem.Chomog/totVol;
        end

        function updateMonitoring(obj,density,step)
            if ~isempty(obj.monitor) == true
                obj.monitor.update(step,{density.fValues});
                obj.monitor.refresh();
                exportgraphics(gcf,'crack.gif','Append',true);
            end
        end

    end
    
end