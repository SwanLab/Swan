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
                paramHole{i} = linspace(0.1,obj.maxParam(i),obj.nSteps(i));
            end
        end
        
        function maxV = computeMaxHoleParams(obj)
            switch obj.holeType
                case 'Circle'
                    maxV = 0.49;
                case 'Square'
                    maxV = 0.98;
                case 'Ellipse'
                    maxV = [0.98,0.98];
                case 'Rectangle'
                    maxV = [0.98,0.98];
                case 'SmoothHexagon'
                    maxV = 0.98;
            end
        end

        function matHomog = computeHomogenization(obj,l)
            mesh = obj.createMesh(l);
            mesh = obj.baseMesh;
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
            switch obj.meshType
                case 'Square'
                    % Dirichlet
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

                    % Periodic
                    sPer{1}.leader = @(coor) isLeft(coor);
                    sPer{1}.follower = @(coor) isRight(coor);
                    sPer{2}.leader = @(coor) isBottom(coor);
                    sPer{2}.follower = @(coor) isTop(coor);
                    sPer{3}.vertex = @(coor) isVertex(coor);
                case 'Hexagon'
                    % Dirichlet
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

                    % Periodic
                    sPer{1}.leader = @(coor) isBottom(coor);
                    sPer{1}.follower = @(coor) isTop(coor);
                    sPer{2}.leader = @(coor) isRightBottom(coor);
                    sPer{2}.follower = @(coor) isLeftTop(coor);
                    sPer{3}.leader = @(coor) isRightTop(coor);
                    sPer{3}.follower = @(coor) isLeftBottom(coor);
                    sPer{4}.vertex = @(coor) isVertex(coor);
            end
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.pointloadFun = [];
            s.periodicFun  = 1;
            s.mesh = mesh;
            bc = BoundaryConditions(s);
            bc.updatePeriodicConditions(sPer);
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
                            apothem   = sqrt(l^2 + (l/2)^2);
                            phi = perimeter*apothem/2;
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