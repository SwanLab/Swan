classdef CoarsePlotSolution < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mesh
        u
        originalBCApplier
        outputFileName
        r
        centroids
    end

    methods
        function obj = CoarsePlotSolution(u, mesh, bcApplier,outputFileName, r, centroids)
            
            if isempty(r)
                obj.plot(u, mesh, bcApplier,outputFileName);
            else
                obj.init(u, mesh, bcApplier,outputFileName,r, centroids);
                obj.makeHole(bcApplier,outputFileName);
                
            end
        end

        % function makeHole(~, x, mesh, bcApplier,outputFileName,r,holeCoords)
        %     x = obj.createXHole(x, mesh, r, holeCoords);
        %     mesh = createHole(mesh, r, holeCoords);
        %     obj.plot(x, mesh, bcApplier,outputFileName);
        % 
        % end

        

    end

    methods (Access = private) %%% treure Static

        function makeHole(obj, bcApplier,outputFileName)
                    
                    %uH = obj.createUHole();
                    meshH = obj.createHole();
                    obj.plot(uH, meshH, bcApplier, outputFileName);
        
        end

        function init(obj, u, mesh, bcApplier,outputFileName,r,centroids)
            obj.mesh        = mesh;
            obj.u           = u;
            obj.originalBCApplier   = bcApplier;
            obj.outputFileName      = outputFileName;
            obj.r                   = r;
            obj.centroids   = centroids;

        end

        function U = createUHole(obj)
            U = obj.u;

            for i = obj.mesh.nnodes:-1:1
                if sqrt( (obj.mesh.coord(i,1) )^2 + (obj.mesh.coord(i,2))^2 ) <= obj.r
                    U(2*i-1:2*i) = [];
                end

            end

        end

        function meshH = createHole(obj)
            % mR    = obj.createReferenceMesh();

            gPar.type         = 'Circle';

            % Initialize the sum as the zero function
            f_sum = @(x) 0;
            u = obj.u;
            mesh = obj.mesh;
            % Loop through each handle and build the sum
            for i = 1:length(obj.r)
                radius = obj.r(i);
                x0 = obj.centroids(i,1);
                y0 = obj.centroids(i,2);
                % f_prev = f_sum;
                f = @(x) (x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2-radius^2;
                s.fHandle = f;
                s.ndimf   = 1;
                s.mesh    = mesh;
                aFun      = AnalyticalFunction(s);
                ls        = aFun.project('P1');

                lsCircle          = ls.fValues;
                lsCircleInclusion = -lsCircle;

                sUm.backgroundMesh = obj.mesh;
                sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
                uMesh              = UnfittedMesh(sUm);
                uMesh.compute(lsCircleInclusion);
                uMeshFun = uMesh.obtainFunctionAtUnfittedMesh(u);
                u = uMeshFun.innerMeshFunction;
                mesh = uMesh.innerMesh;
                %f = @(x) (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)<radius) + ...
                %        (sqrt((x(1,:,:)-x0).^2+(x(2,:,:)-y0).^2)>=radius) ; % Save previous sum to avoid recursion
                % f_sum = @(x) f_prev(x) + f(x);  % Add new handle
            end

            s.fHandle = f_sum;
            s.ndimf   = 1;
            s.mesh    = obj.mesh;
            aFun      = AnalyticalFunction(s);
            ls        = aFun.project('P1');

            for i = 1:size(obj.r,2)
                gPar.radius       = obj.r(i);
                gPar.xCoorCenter  = obj.centroids(i,1);
                gPar.yCoorCenter  = obj.centroids(i,2);
                g                 = GeometricalFunction(gPar);
            end
            lvSet = obj.createLevelSetFunction(mR);
            uMesh = obj.computeUnfittedMesh(mR,lvSet);
            meshH  = uMesh.createInnerMesh();


        end

        function newMesh = createReferenceMesh(obj)
            xmin = -1;
            xmax = 1;
            ymin = -1;
            ymax = 1;

            n1 = numel(find(obj.mesh.coord(:,1)==xmin));
            n2 = n1;%numel(find(obj.mesh.coord(:,2)==ymin));

           
             %UnitMesh better
            x1      = linspace(xmin,xmax,n1);
            x2      = linspace(ymin,ymax,n2);
            [xv,yv] = meshgrid(x1,x2);
            [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
            s.coord  = V(:,1:2);
            s.connec = F;
            newMesh = Mesh.create(s);
        end

        function levelSet = createLevelSetFunction(obj,mR)
            sLS.type        = 'CircleInclusion';
            sLS.xCoorCenter = 0;
            sLS.yCoorCenter = 0;
            sLS.radius      = obj.r;
            g               = GeometricalFunction(sLS);
            lsFun           = g.computeLevelSetFunction(mR);
            levelSet        = lsFun.fValues;
        end

        function uMesh = computeUnfittedMesh(obj, mR,levelSet)
            sUm.backgroundMesh = mR;
            sUm.boundaryMesh   = mR.createBoundaryMesh();
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(levelSet);
        end

        function plot(~, x, mesh, bcApplier,outputFileName)
            row = 0;
            col = 0;
            iter = 0;
            flag = 0;

            if ~isempty(bcApplier)
                x = bcApplier.reducedToFullVectorDirichlet(x);
            end
            if nargin <7
                flag =0;
            end
            %             xFull = bc.reducedToFullVector(x);
            if size(x,2)==1
                s.fValues = reshape(x,2,[])';
            else
                s.fValues = x;
            end
            %

            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 2;
            s.order = 'P1';
            xF = LagrangianFunction(s);
            %             xF.plot();
            if flag == 0
                xF.print(['domain',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 1
                xF.print(['DomainResidual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 2
                xF.print(['Residual',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 3
                xF.print(['domainFine',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            elseif flag == 4
                xF.print(['domainNeuman',num2str(row),num2str(col),'_',num2str(iter)],'Paraview')
            end
                        
            fileName = ['domain',num2str(row),num2str(col),'_',num2str(iter)];

            s = dir(pwd);
            s = struct2table(s);
            idx = startsWith(s.name, fileName);
            s = s(idx,:);
            oldFileName = s.name;
            newFileName = replace(oldFileName, fileName, outputFileName);
            fclose('all');
            movefile(oldFileName{1}, newFileName{1})

        end
        

    end
end