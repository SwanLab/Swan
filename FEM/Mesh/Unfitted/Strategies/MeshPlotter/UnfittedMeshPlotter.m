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
            hold on
            obj.plotAll();            
        end
        
        function plotBoundary(obj)
            figure
            hold on
            obj.plotMesh(obj.uMesh.backgroundMesh);
            obj.plotMesh(obj.uMesh.boundaryCutMesh);
            obj.plotUnfittedBoundaryMesh();
        end  
        
        function plotAll(obj)
            uM = obj.uMesh;
            obj.plotMesh(uM.backgroundMesh);
            obj.plotMesh(uM.innerMesh);
            obj.plotMesh(uM.innerCutMesh);
            obj.plotMesh(uM.boundaryCutMesh);            
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.uMesh = cParams.uMesh;
        end
        
        function plotUnfittedBoundaryMesh(obj)  
            uB = obj.uMesh.unfittedBoundaryMesh;
            uBoundaryMeshes = uB.getActiveMesh();
            for imesh = 1:numel(uBoundaryMeshes)
                uM = uBoundaryMeshes{imesh};
                uM.plotAll();
            end            
        end       
        
    end
    
    methods (Access = private, Static)
        
        function plotMesh(mesh)
            s.mesh = mesh;
            mP = MeshPlotter(s);
            mP.plot();
        end
        
    end    
    
end