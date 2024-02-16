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
            obj.initVol(cParams)
        end
        
        function hMin = computeMinCellSize(obj)
            % wrong
            x1(:,1) = obj.coord(obj.connec(:,1),1);
            x1(:,2) = obj.coord(obj.connec(:,1),2);
            x2(:,1) = obj.coord(obj.connec(:,2),1);
            x2(:,2) = obj.coord(obj.connec(:,2),2);
            x3(:,1) = obj.coord(obj.connec(:,3),1);
            x3(:,2) = obj.coord(obj.connec(:,3),2);
            x1x2 = (x2-x1);
            x2x3 = (x3-x2);
            x1x3 = (x1-x3);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
            n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
            hs = min([n12,n23,n13],[],2);
            hMin = min(hs);
        end

        function hMean = computeMeanCellSize(obj)
            % wrong
            x1(:,1) = obj.coord(obj.connec(:,1),1);
            x1(:,2) = obj.coord(obj.connec(:,1),2);
            x2(:,1) = obj.coord(obj.connec(:,2),1);
            x2(:,2) = obj.coord(obj.connec(:,2),2);
            x3(:,1) = obj.coord(obj.connec(:,3),1);
            x3(:,2) = obj.coord(obj.connec(:,3),2);
            x1x2 = (x2-x1);
            x2x3 = (x3-x2);
            x1x3 = (x1-x3);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
            n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
            hs = max([n12,n23,n13],[],2);
            hMean = max(hs);
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
    
    methods (Access = private)
        
        function initVol(obj,cParams)
        end
        
    end
    
end