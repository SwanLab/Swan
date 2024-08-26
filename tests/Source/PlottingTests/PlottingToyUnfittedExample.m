classdef PlottingToyUnfittedExample < handle
    
    properties (SetAccess = protected, GetAccess = private)
        coord
        connec
        levelSet
        testName
    end
    
    properties (Access = private)
       backgroundMesh 
       boundaryMeshes
       unfittedMesh
    end
    
    methods (Access = public)
        
        function obj = PlottingToyUnfittedExample(cParams)
            figure()
            obj.init(cParams);
            obj.compute();
            obj.plot();
        end

        function passed = hasPassed(obj)
            d = load(obj.testName);
            unfittedMesh = obj.unfittedMesh;
            itIs = isequaln(unfittedMesh,d.unfittedMesh);
            passed = itIs;
%             save(obj.testName, 'unfittedMesh', '-append')
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.coord    = cParams.coord;
            obj.connec   = cParams.connec;
            obj.levelSet = cParams.levelSet;
            obj.testName = cParams.testName;
        end

        function compute(obj)
            obj.computeBackgroundMesh();
            obj.computeBoundaryMeshes();
            obj.computeUnfittedMesh();
        end
        
        function plot(obj)
            figure();
            obj.plotUnfittedMesh();
            obj.plotGaussPointsInUnfittedMesh(obj.unfittedMesh)
            obj.plotUnfittedBoundaryMesh();
        end

        function computeBackgroundMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            obj.backgroundMesh = Mesh.create(s);
        end
        
        function computeBoundaryMeshes(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.dimension = 1:obj.backgroundMesh.ndim;
            bC = BoundaryMeshCreatorFromRectangularBox(s);
            obj.boundaryMeshes = bC.create();
        end
        
        function computeUnfittedMesh(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.boundaryMesh   = obj.boundaryMeshes;
            uM = UnfittedMesh(s);
            uM.compute(obj.levelSet);
            obj.unfittedMesh = uM;
        end
        
        function plotUnfittedMesh(obj)
            obj.unfittedMesh.plot();
        end
        
        function plotGaussPointsInUnfittedMesh(obj,uMesh)
            obj.plotInnerGaussPoints(uMesh.innerMesh);
            obj.plotInnerCutGaussPoints(uMesh.innerCutMesh);
            obj.plotBoundaryCutGaussPoints(uMesh.boundaryCutMesh);
        end
        
        function plotInnerGaussPoints(obj,innerMesh)
             if ~isempty(innerMesh)
                mesh = innerMesh.mesh;
                color = [0.8500    0.3250    0.0980];
                obj.plotGaussPoints(mesh,2,color)
             end
        end
        
        function plotInnerCutGaussPoints(obj,innerCutMesh)
            if ~isempty(innerCutMesh)
                mesh = innerCutMesh.mesh;
                color = [0.4940    0.1840    0.5560];
                obj.plotGaussPoints(mesh,2,color)
            end
        end
        
        function plotBoundaryCutGaussPoints(obj,boundaryCutMesh)
            if ~isempty(boundaryCutMesh)
                mesh = boundaryCutMesh.mesh;
                color = [0.9290    0.6940    0.1250];
                obj.plotGaussPoints(mesh,2,color);
            end            
        end
        
        function plotUnfittedBoundaryMesh(obj)
            bMesh = obj.unfittedMesh.unfittedBoundaryMesh;
            meshes = bMesh.getActiveMesh();
            for iMesh = 1:numel(meshes)
                mesh = meshes{iMesh};
                obj.plotGaussPointsInUnfittedMesh(mesh)
            end
        end

    end

    methods (Access = private, Static)

        function plotGaussPoints(mesh,quadType,color)
            quad = Quadrature.create(mesh,quadType);
            xCutB = mesh.computeXgauss(quad.posgp);
            p = plot(xCutB(1,:),xCutB(2,:),'s','MarkerSize',5);
            p.MarkerEdgeColor = color;
            p.Color = color;
            p.MarkerFaceColor = color;
        end

    end
end