classdef ElasticProblemPunzonV2 < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
        filename
    end

    methods (Access = public)

        function obj = ElasticProblemPunzonV2()
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

        % function computeElasticProperties(obj)
        %     E  = 1;
        %     nu = 1/3;
        %     obj.young   = ConstantFunction.create(E,obj.mesh);
        %     obj.poisson = ConstantFunction.create(nu,obj.mesh);
        % end

        % function createMaterial(obj)
        %     s.type    = 'ISOTROPIC';
        %     s.ptype   = 'ELASTIC';
        %     s.ndim    = obj.mesh.ndim;
        %     s.young   = obj.young;
        %     s.poisson = obj.poisson;
        %     tensor    = Material.create(s);
        %     obj.material = tensor;
        % end

        % function createElasticProblem(obj)
        %     file = 'punzon';
        %     a.fileName = file;
        %     s = FemDataContainer(a)
        %     s.mesh = obj.mesh;
        %     s.scale = 'MACRO';
        %     s.material = obj.material;
        %     s.dim = '2D';
        %     s.boundaryConditions = obj.createBoundaryConditions();
        %     s.solverType = 'REDUCED';
        %     s.solverMode = 'DISP';
        %     s.solverCase = DirectSolver();
        %     fem = ElasticProblem(s);
        %     fem.solve();
        %     fem.uFun.print('results_fem_dispPunzon', 'Paraview') % print using Paraview
        %     obj.stateProblem = fem;
        % end

        function bc = createBoundaryConditions(obj)
            femReader = FemInputReaderGiD();
            s         = femReader.read(obj.filename);
            sPL       = obj.computeCondition(s.pointload);
            sDir      = obj.computeCondition(s.dirichlet);

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            pointloadFun = [];
            for i = 1:numel(sPL)
                pl = PointLoad(obj.mesh, sPL{i});
                pointloadFun = [pointloadFun, pl];
            end
            s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh         = obj.mesh;
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