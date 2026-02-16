classdef ElasticProblemPunzon < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
        filename
    end

    methods (Access = public)

        function obj = ElasticProblemPunzon()
            obj.init();
            a.fileName = obj.filename;
            s = FemDataContainer(a);
            obj.createMesh(s);
            s.boundaryConditions=obj.createBoundaryConditions();
            % obj.computeElasticProperties();
            % obj.createMaterial();
            %obj.createElasticProblem();
            fem = PhysicalProblem.create(s);
            fem.solve();
            fem.print('prova1','Paraview');
        end

    end

    methods (Access = private)

        function init(obj)
            obj.filename='punzon';
        end

        function createMesh(obj,s)
            m = s.mesh;
            con = m.connec;
            q = Quadrature.create(m,0);
            dv = m.computeDvolume(q);
            negElem=find(dv<=0);
            con(negElem,:) = [];
            sM.coord = m.coord;
            sM.connec = con;
            m2 = Mesh.create(sM);
            m2 = m2.computeCanonicalMesh();
            obj.mesh = m2;
        end



        function bc = createBoundaryConditions(obj)
            xMin = min(obj.mesh.coord(:,1));
            xMax = max(obj.mesh.coord(:,1));
            yMin = min(obj.mesh.coord(:,2));
            yMax = max(obj.mesh.coord(:,2));
            zMin = min(obj.mesh.coord(:,3));
            zMax = max(obj.mesh.coord(:,3));

            isBottom = @(coor) abs(coor(:,3) - zMin) < 1e-6; % cara llisa
            isTop = @(coor) abs(coor(:,3) - zMax) < 1e-6; % tornillos
            % isGuide1 = @(coor) abs(coor(:,1) - xMin) < 1e-6 & abs(coor(:,2) - yMin) < 1e-6;
            % isGuide2 = @(coor) abs(coor(:,1) - xMin) < 1e-6 & abs(coor(:,2) - yMax) < 1e-6;
            % isGuide3 = @(coor) abs(coor(:,1) - xMax) < 1e-6 & abs(coor(:,2) - yMin) < 1e-6;
            % isGuide4 = @(coor) abs(coor(:,1) - xMax) < 1e-6 & abs(coor(:,2) - yMax) < 1e-6;

            isForce = @(x) (x(3,:,:) - zMin) < 1e-2;

            sDir{1}.domain    = @(coor) isTop(coor);
            sDir{1}.direction = [1,2,3];
            sDir{1}.value     = 0;

            sPL{1}.domain    = isForce;
            [bMesh,~]=obj.mesh.createSingleBoundaryMesh();
            sPL{1}.fun       = ConstantFunction.create([0,0,1],bMesh);


            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = TractionLoad(obj.mesh, sPL{i}, 'FUNCTION');
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun = [];
            s.mesh        = obj.mesh;

            bc = BoundaryConditions(s);
    
        end

    end

    methods (Static, Access=private)
        function sCond = computeCondition(conditions)
            nodes = @(coor) 1:size(coor,1);
            dirs  = unique(conditions(:,2));
            j     = 0;
            for k = 1:length(dirs)
                rowsDirk = ismember(conditions(:,2),dirs(k));
                u        = unique(conditions(rowsDirk,3));
                for i = 1:length(u)
                    rows   = conditions(:,3)==u(i) & rowsDirk;
                    isCond = @(coor) ismember(nodes(coor),conditions(rows,1));
                    j      = j+1;
                    sCond{j}.domain    = @(coor) isCond(coor);
                    sCond{j}.direction = dirs(k);
                    sCond{j}.value     = u(i);
                end
            end
        end
    end

end