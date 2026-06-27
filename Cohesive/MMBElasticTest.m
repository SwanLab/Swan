classdef MMBElasticTest < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
        leverDisp
    end

    methods (Access = public)

        function obj = MMBElasticTest()
            obj.init();
            obj.createMesh();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
            obj.postProcess();
        end

    end

    methods (Access = private)

        function init(obj)
            obj.leverDisp = 0.6;
        end

        function createMesh(obj)
            a.fileName = 'CThickLeverTriangular';
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function computeElasticProperties(obj)
            E  = 120e3;
            nu = 0.3;
            obj.young   = ConstantFunction.create(E,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = obj.young;
            s.poisson = obj.poisson;
            tensor    = Material.create(s);
            obj.material = tensor;
        end

        function solveElasticProblem(obj)
            s.mesh = obj.mesh;
            s.scale = 'MACRO';
            s.material = obj.material;
            s.dim = '2D';
            s.boundaryConditions = obj.createBoundaryConditions();
            s.solverType = 'REDUCED';
            s.solverMode = 'DISP';
            s.solverCase = 'DIRECT';
            fem = ElasticProblem(s);
            fem.solve();
            obj.stateProblem = fem;
        end

        function bc = createBoundaryConditions(obj)
            isUp   = @(coor) abs(coor(:,2) - max(coor(:,2))) < 1e-12;
            isDown = @(coor) abs(coor(:,2) - min(coor(:,2))) < 1e-12;
            isLeft = @(coor)  abs(coor(:,1)-min(coor(:,1))) < 1e-12;
            isRight = @(coor)  abs(coor(:,1)-max(coor(:,1))) < 1e-12;
           
            sDir.domain    = @(coor) isDown(coor) & isLeft(coor);
            sDir.direction = [1,2];
            sDir.value     = 0;
            Dir1 = DirichletCondition(obj.mesh,sDir);

            sDir.domain    = @(coor) isDown(coor) & isRight(coor);
            sDir.direction = [2];
            sDir.value     = 0;
            Dir2 = DirichletCondition(obj.mesh,sDir);

            isLeftLever    = @(coor) abs(coor(:,1) - 11.82) < 1e-12;
            sDir.domain    = @(coor) isUp(coor) & isLeftLever(coor);
            sDir.direction = [2];
            sDir.value     = -obj.leverDisp;
            Dir3 = DirichletCondition(obj.mesh,sDir);

            s.mesh = obj.mesh;
            s.dirichletFun = [Dir1 Dir2 Dir3];
            s.pointloadFun = [];
            s.periodicFun = [];
            bc = BoundaryConditions(s);
        end

        function postProcess(obj)
            uFun = obj.stateProblem.uFun;

            forces = obj.stateProblem.forces;
            isLeftLever      =  abs(obj.mesh.coord(:,1) - 11.82) < 1e-12;
            isUp   = abs(obj.mesh.coord(:,2) - max(obj.mesh.coord(:,2))) < 1e-12;
            isLeverHandle    =  isUp & isLeftLever;
            nodes = find(isLeverHandle);
            dofsY= (nodes-1)*uFun.ndimf + 2;
            totReact = abs(sum(forces(dofsY)));

            disp(['Total Reaction Force in Y: ', num2str(totReact)]);
            disp(['Lever Displacement: ', num2str(obj.leverDisp)]);

            obj.writeGiDPostRes(uFun);
        end


        function writeGiDPostRes(obj,uFun)
            fname = fullfile('CThickLeverTriangular.post.res');
            fid = fopen(fname,'w');
            fprintf(fid,'GiD Post Results File 1.0\n');
            fprintf(fid,'Result "Displacement" "Static" 1 Vector OnNodes\n');
            fprintf(fid,'ComponentNames "Ux" "Uy"\n');
            fprintf(fid,'Values\n');
            for i = 1:size(obj.mesh.coord,1)
                fprintf(fid,'%d %.16e %.16e\n',i,uFun.fValues(i,1),uFun.fValues(i,2));
            end
            fprintf(fid,'End Values\n');
            fclose(fid);
        end

    end

end