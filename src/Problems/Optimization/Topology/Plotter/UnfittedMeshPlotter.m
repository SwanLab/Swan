classdef UnfittedMeshPlotter < handle
    
    properties (Access = private)
        uMesh
        faceColor
        edgeAlpha
        faceAlpha
    end
    
    methods (Access = public)
        
        function obj = UnfittedMeshPlotter(cParams)
            obj.init(cParams)
        end
        
        function plotDomainInColor(obj,color)
            obj.faceColor = color;
            obj.edgeAlpha = 1;
            obj.faceAlpha = 0.9;
            obj.plotInner();
            obj.plotInnerCut();
        end
        
        function plotDomain(obj)
            obj.faceColor = 'red';
            obj.edgeAlpha = 0.5;
            obj.faceAlpha = 0.3;
            obj.plotAll();
            obj.addLighting();
        end
      
        function plotBoundary(obj)
            obj.plotBackground();
            obj.plotBoundaryCutMesh();
            obj.plotUnfittedBoundaryMesh();
            obj.addLighting();
        end
        
        function plotAll(obj)
            obj.plotBackground();
            obj.plotInner();
            obj.plotInnerCut();
            obj.plotBoundaryCutMesh();
            obj.plotUnfittedBoundaryMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.uMesh = cParams.uMesh;
        end
        
        function plotBackground(obj)
            s.mesh = obj.uMesh.backgroundMesh;
            s.isBackground = true;
            s.faceColor = 'red';
            s.faceAlpha = 0.3;
            s.edgeAlpha = 0.5;
            obj.plotMesh(s);
        end
        
        function plotInner(obj)
            m = obj.uMesh.innerMesh;
            obj.plotSubMesh(m);
        end
        
        function plotInnerCut(obj)
            m = obj.uMesh.innerCutMesh;
            obj.plotSubMesh(m);
        end
        
        function plotBoundaryCutMesh(obj)
            m = obj.uMesh.boundaryCutMesh;
            obj.plotSubMesh(m);
        end
        
        function plotUnfittedBoundaryMesh(obj)
            uB = obj.uMesh.unfittedBoundaryMesh;
            if ~isempty(uB)
                uBoundaryMeshes = uB.getActiveMesh();
                for imesh = 1:numel(uBoundaryMeshes)
                    uM = uBoundaryMeshes{imesh};
                    uM.plotAll();
                end
            end
        end
        
        function plotSubMesh(obj,uM)
            s.faceColor = obj.faceColor;
            s.edgeAlpha = obj.edgeAlpha;
            s.faceAlpha = obj.faceAlpha;
            s.isBackground = false;
            if ~isempty(uM)
                s.mesh = uM.mesh;
                obj.plotMesh(s);
            end
        end
        
        
        function addLighting(obj)
            delete(findall(gcf,'Type','light'))
            if ~isempty(obj.uMesh.innerMesh)
                if obj.uMesh.innerMesh.mesh.ndim == 3
                    l = lightangle(120,30);
                    material metal
                    l = lightangle(240,30);
                    material metal
                    l = lightangle(90,180);
                    material metal
                    l = lightangle(90,0);
                    material metal
                    view([1 1 1])
                end
            end
        end
        
       
        
    end
    
    
    methods (Access = private, Static)
        
        function plotMesh(s)
            if ~isempty(s.mesh)
                mP = MeshPlotter(s);
                mP.plot();
            end
        end
    end
    
end