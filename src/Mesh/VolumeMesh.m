classdef VolumeMesh < Mesh
    
    properties (Access = public)
        geometryType = 'Volume';
    end
    
    properties (Access = private)
        cParams;
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = VolumeMesh(cParams)
            obj = obj@Mesh(cParams);
        end

        function plot(obj)
            gPar.type         = 'Full';
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(obj);
            lsCircle          = phiFun.fValues;
            
            sUm.backgroundMesh = obj;
            sUm.boundaryMesh   = obj.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsCircle);
            uMesh.plot
        end
        
    end
       
end