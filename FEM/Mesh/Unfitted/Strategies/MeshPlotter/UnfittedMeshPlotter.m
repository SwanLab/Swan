classdef UnfittedMeshPlotter < handle
    
    properties (Access = private)
        uMesh
    end
    
    methods (Access = public)
        
        function obj = UnfittedMeshPlotter(cParams)
            obj.init(cParams)
        end
        
        function plotDomain(obj)
            figure
            clf            
            hold on
            obj.plotAll();
        end
        
        function plotBoundary(obj)
            figure(1)
            clf
            hold on
            obj.plotBackground();
            obj.plotBoundaryCutMesh();
            obj.plotUnfittedBoundaryMesh();
        end
        
        function plotAll(obj)
            obj.plotBackground();
            obj.plotInner();
            obj.plotInnerCut();
            obj.plotBoundaryCutMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.uMesh = cParams.uMesh;
        end
        
        function plotBackground(obj)
            s.mesh = obj.uMesh.backgroundMesh;
            s.isBackground = true;
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
            uBoundaryMeshes = uB.getActiveMesh();
            for imesh = 1:numel(uBoundaryMeshes)
                uM = uBoundaryMeshes{imesh};
                uM.plotAll();
            end
        end
        
        function plotSubMesh(obj,uM)
            if ~isempty(uM)
                s.mesh = uM.mesh;
                obj.plotMesh(s);
            end
        end        
        
    end
    
    methods (Access = private, Static)
        
        function plotMesh(s)
            s = SettingsMeshPlotter(s);
            if ~isempty(s.mesh)
                mP = MeshPlotter(s);
                mP.plot();
            end
        end
    end
    
end