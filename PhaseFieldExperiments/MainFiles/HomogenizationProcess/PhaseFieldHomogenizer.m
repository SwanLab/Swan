classdef PhaseFieldHomogenizer < handle

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
        maxParam
    end

    methods (Access = public)
        
        function obj = PhaseFieldHomogenizer(cParams)
            obj.init(cParams);
            obj.defineMesh();
        end
        
        function [mat,phi] = computeHomogMaterial(obj)
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
            mat = obj.assembleResults(mat);
            phi = obj.assembleResults(phi);
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
            s.coord = MC.coord;
            s.connec = MC.connec;
            obj.baseMesh = Mesh.create(s);
        end

        function paramHole = computeHoleParams(obj)
            obj.maxParam = obj.computeMaxHoleParams();
            nParam = length(obj.maxParam);
            paramHole = cell(1,nParam);
            for i=1:nParam
                paramHole{i} = linspace(0.01,obj.maxParam(i),obj.nSteps(i));
            end
        end
        
        function maxV = computeMaxHoleParams(obj)
            switch obj.holeType
                case 'Circle'
                    maxV = 0.5;
                case 'Ellipse'
                    maxV = [0.5,0.5];
                case 'Square'
                    maxV = 0.98;
                case 'Rectangle'
                    maxV = [0.98,0.98];
                case 'SmoothHexagon'
                    maxV = [0.98];
            end
        end

        function matHomog = computeHomogenization(obj,l)
            mesh = obj.createMesh(l);
            mat = obj.createMaterial(mesh);
            matHomog = obj.solveElasticMicroProblem(mesh,mat);
        end

        function mesh = createMesh(obj,l)
            ls = obj.computeLevelSet(obj.baseMesh,l);
            sUm.backgroundMesh = obj.baseMesh;
            sUm.boundaryMesh   = obj.baseMesh.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(ls);
            holeMesh = uMesh.createInnerMesh();
            figure()
            holeMesh.plot
            mesh = holeMesh;
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
                case 'Ellipse'
                    gPar.xSide  = l(1);
                    gPar.ySide  = l(2);
                case 'Square'
                    gPar.radius = l;
                case 'Rectangle'
                    gPar.xSide  = l(1);
                    gPar.ySide  = l(2);
                case 'Hexagon'
                    gPar.radius = l;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
                case 'SmoothHexagon'
                    gPar.radius = l;
                    gPar.normal = [0 1; sqrt(3)/2 1/2; sqrt(3)/2 -1/2];
            end
            g                  = GeometricalFunction(gPar);
            phiFun             = g.computeLevelSetFunction(mesh);
            lsCircle           = phiFun.fValues;
            ls = -lsCircle;
        end

        function mat = createMaterial(obj,mesh)
            young   = ConstantFunction.create(obj.E,mesh);
            poisson = ConstantFunction.create(obj.nu,mesh);
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.mesh    = mesh;
            s.young   = young;
            s.poisson = poisson;
            mat    = Material.create(s);
        end

        function matHomog = solveElasticMicroProblem(obj,mesh,material)
            figure(1)
            cla reset
            mesh.plot

            s.mesh = mesh;
            s.material = material;
            s.scale = 'MICRO';
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions(mesh);
            s.solverCase = 'DIRECT';
            s.solverType = 'REDUCED';
            s.solverMode = 'FLUC';
            fem = ElasticProblemMicro(s);
            fem.solve();
            matHomog = fem.Chomog;
        end

        function bc = createBoundaryConditions(obj,mesh)
            isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
            isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
            isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);
            isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
            
            % Dirichlet
            
            sDir{1}.domain    = @(coor) isTop(coor) & isLeft(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;
            
            sDir{2}.domain    = @(coor) isTop(coor) & isRight(coor);
            sDir{2}.direction = [1,2];
            sDir{2}.value     = 0;
            
            sDir{3}.domain    = @(coor) isBottom(coor) & isLeft(coor);
            sDir{3}.direction = [1,2];
            sDir{3}.value     = 0;
            
            sDir{4}.domain    = @(coor) isBottom(coor) & isRight(coor);
            sDir{4}.direction = [1,2];
            sDir{4}.value     = 0;
            
            % Periodic (NOTE: actually sets all boundaries as periodic hehe)
            
            sPer{1}.leader = @(coor) isLeft(coor);
            sPer{1}.follower = @(coor) isRight(coor);

            dirichletFun = [];
            periodicFun = [];

            for i = 1:numel(sDir)
                dir = DirichletCondition(mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end

            for i = 1:numel(sPer)
                per = PeriodicCondition(mesh, sPer{i});
                periodicFun = [periodicFun, per];
            end

            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = periodicFun;
            s.mesh = mesh;
            bc = BoundaryConditions(s);
        end
        
        function phi = computeDamageMetric(obj,l)
            max = obj.maxParam;
            switch obj.damageType
                case 'Area'
                    switch obj.holeType
                        case 'Circle'
                            phi = l^2/(max^2);
                        case 'Ellipse'
                            phi = (l(1)*l(2))/(max(1)*max(2));
                        case 'Square'
                            phi = l^2/max^2;
                        case 'Rectangle'
                            phi = (l(1)*l(2))/(max(1)*max(2));
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