classdef FittedNaca < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        dY
        stateProblem
        filename
        nacaNodes
        mNew
    end

    methods (Access = public)

        function obj = FittedNaca()
            obj.init();
            obj.createMesh();
            obj.readNacaNodes();
            obj.computeUpdatedNaca();
            obj.computeElasticProperties();
            obj.createMaterial();
            obj.solveElasticProblem();
            obj.createNewMesh();
        end

    end

    methods (Access = private)

        function init(obj)
            file = 'Naca';
            obj.filename = file;
        end

        function createMesh(obj)
            a.fileName = obj.filename;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end

        function readNacaNodes(obj)
            a.fileName = obj.filename;
            s = FemDataContainer(a);            
            obj.nacaNodes = unique(s.bc.dirichlet(:,1));
        end        

        function computeElasticProperties(obj)
            E  = 1;
            nu = 1/3;
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

        function computeUpdatedNaca(obj)
            mRef      = 0.02;
            pRef      = 0.4;
            tRef      = 0.12;
            isUP = obj.compueIsUp(mRef,pRef,tRef);

            m      = 0.01;
            p      = 0.2;
            t      = 0.05;
            xNaca = obj.mesh.coord(obj.nacaNodes,1);
            yNaca = obj.mesh.coord(obj.nacaNodes,2);                                    
            [yu,yl,~] = obj.computeNacaFunction(m,p,t,xNaca);
                
            incY(isUP)  = yu(isUP) - yNaca(isUP);
            incY(~isUP) = -yl(~isUP) + yNaca(~isUP);
            
            obj.dY = incY;
        end

        function isUP = compueIsUp(obj,m,p,t)              
            xNaca = obj.mesh.coord(obj.nacaNodes,1);
            yNaca = obj.mesh.coord(obj.nacaNodes,2);                        
            [~,~,yc] = obj.computeNacaFunction(m,p,t,xNaca);
            isUP = yNaca > yc;
        end

        function bc = createBoundaryConditions(obj)
            dYCond      = obj.nacaNodes;
            dYCond(:,2) = 2;
            dYCond(:,3) = obj.dY;
            sDir        = obj.computeCondition(dYCond);

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            s.pointloadFun = [];
            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

        function createNewMesh(obj)
            u = obj.stateProblem.uFun;

            s.coord  = obj.mesh.coord + u.fValues;
            s.connec = obj.mesh.connec;
            obj.mNew = Mesh.create(s);
        end

    end

    methods (Static, Access=private)
        function [yu,yl,yc] = computeNacaFunction(m,p,t,xNaca)
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));

            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            yl    = yc - yt.*cos(theta);

        end

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