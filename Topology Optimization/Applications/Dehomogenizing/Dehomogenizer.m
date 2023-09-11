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
                figure(1)
                clf
                uM = obj.createUnfittedMesh(ls{i});
                uM.plotStructureInColor('black');
                drawnow
                %    uM.plotComponents();
             %   obj.saveImage()
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
            L = obj.mesh.computeCharacteristicLength();
            obj.epsilons = L./obj.nCells;
            obj.epsilons
        end            

        function o = computeOrientedMappingComputer(obj)
            s.orientationP0 = obj.theta;
            s.mesh  = obj.mesh;
            o = OrientedMappingComputer(s);
        end

        function computeLevelSet(obj)
            s.type               = 'periodicAndOriented';
            s.mesh               = obj.mesh;
            s.orientationVectors = obj.computeOrientedMappingComputer();
            s.cellLevelSetParams = obj.cellLevelSetParams;
            lSet = LevelSetCreator.create(s);
            ls = lSet.computeLS(obj.epsilons);
            obj.levelSet = ls;  
            obj.fineMesh = lSet.getFineMesh();
        end

        function uM = createUnfittedMesh(obj,ls)
            s.boundaryMesh   = obj.createBoundaryMesh();
            s.backgroundMesh = obj.fineMesh;
            uM = UnfittedMesh(s);
            uM.compute(ls);
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