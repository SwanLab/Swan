classdef Dehomogenizer < handle

    properties (Access = private)
        levelSet
        fineMesh
        epsilons
    end

    properties (Access = private)
        cellLevelSetParams
        nCells
        theta
        mesh
        remesher
    end

    methods (Access = public)

        function obj = Dehomogenizer(cParams)
            obj.init(cParams);
        end

        function ls = compute(obj)
            obj.createEpsilons();
            obj.computeLevelSet();
            ls = obj.levelSet;            
        end

        function plot(obj)
            ls = obj.levelSet;
            for i = 1:numel(ls)
                fh = figure(1);
                clf
                uM = obj.createUnfittedMesh(ls{i});
                uM.plotStructureInColor('black');
                drawnow
                %uM.plotComponents();
             %   obj.saveImage()
           %     fh.WindowState = 'maximized';
           %     axis off
           %    saveas(gcf,['/home/alex/Desktop/CantileverDehomog/Example',num2str(i),'.png'])

            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nCells             = cParams.nCells;
            obj.theta              = cParams.theta;
            obj.cellLevelSetParams = cParams.cellLevelSetParams;
            obj.mesh               = cParams.mesh;
        end

        function createEpsilons(obj)
            xmax = max(obj.mesh.coord(:,1));
            xmin = min(obj.mesh.coord(:,2));
            L    = xmax - xmin;
            e = L./obj.nCells;
            obj.epsilons = (1 + 0.01*linspace(0,100,1))*e/2;
        end            

        function o = computeOrientedMappingComputer(obj)
            s.orientationP1 = obj.theta;
            s.mesh  = obj.mesh;
            o = OrientedMappingComputer(s);
        end

        function computeLevelSet(obj)
       %     s.type               = 'PeriodicAndOriented';
            s.mesh               = obj.mesh;
            s.orientationVectors = obj.computeOrientedMappingComputer();
            s.cellLevelSetParams = obj.cellLevelSetParams;
            ls                   = LevelSetPeriodicAndOriented(s);
            obj.levelSet = ls.computeLS(obj.epsilons);  
            obj.fineMesh = obj.mesh.remesh();%ls.getFineMesh();
        end

        function uM = createUnfittedMesh(obj,ls)
            s.boundaryMesh   = obj.createBoundaryMesh();
            s.backgroundMesh = obj.fineMesh;
            uM = UnfittedMesh(s);
            uM.compute(ls.fValues);
        end

        function bM = createBoundaryMesh(obj)
            sB.backgroundMesh = obj.fineMesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bMc = BoundaryMeshCreator.create(sB);
            bM  = bMc.create();
        end        

        function saveImage(obj)
            xmin = min(obj.mesh.coord(:,1));
            xmax = max(obj.mesh.coord(:,1));
            ymin = min(obj.mesh.coord(:,2));
            ymax = max(obj.mesh.coord(:,2));
            axis([xmin xmax ymin ymax])
            set(gca, 'Visible', 'off')
            exportgraphics(gcf,'testAnimated2.gif','Append',true);
            close all
        end

    end
end