classdef FittedNaca < handle

    properties (Access = private)
        mesh
        young
        poisson
        material
        stateProblem
        filename
        nacaNodes
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

            m      = 0.02;
            p      = 0.4;
            t      = 0.12;
            xNaca = obj.mesh.coord(obj.nacaNodes,1);
            yNaca = obj.mesh.coord(obj.nacaNodes,2);                                    
            [yu,yl,yc] = obj.computeNacaFunction(mRef,pRef,tRef,xNaca);
                
            incY(isUP)  = yu - yNaca(isUP);
            incY(~isUP) = yl + yNaca(isUP);

        end

        function [yu,yl,yc] = computeNacaFunction(obj,m,p,t,xNaca)

            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                (xNaca>p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1036*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                (xNaca>p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));

            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            yl    = yc - yt.*cos(theta);

        end

        function isUP = compueIsUp(obj,m,p,t)              
            xNaca = obj.mesh.coord(obj.nacaNodes,1);
            yNaca = obj.mesh.coord(obj.nacaNodes,2);                        
            [yu,yl,yc] = obj.computeNacaFunction(m,p,t,xNaca);
            isUP = yNaca > yc;
        end

        function bc = createBoundaryConditions(obj)
            xMax    = max(obj.mesh.coord(:,1));
            yMax    = max(obj.mesh.coord(:,2));
            isDir   = @(coor)  abs(coor(:,1))==0;
            isForce = @(coor)  (abs(coor(:,1))==xMax & abs(coor(:,2))>=0.3*yMax & abs(coor(:,2))<=0.7*yMax);

            sDir{1}.domain    = @(coor) isDir(coor);
            sDir{1}.direction = [1,2];
            sDir{1}.value     = 0;

            sPL{1}.domain    = @(coor) isForce(coor);
            sPL{1}.direction = 2;
            sPL{1}.value     = -1;

            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;

            % pointloadFun = [];
            % for i = 1:numel(sPL)
            %     pl = PointLoad(obj.mesh, sPL{i});
            %     pointloadFun = [pointloadFun, pl];
            % end
            % s.pointloadFun = pointloadFun;

            s.periodicFun  = [];
            s.mesh = obj.mesh;
            bc = BoundaryConditions(s);
        end

    end

end