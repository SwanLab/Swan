classdef SurfaceMesh < Mesh
    
    properties (Access = public)
        geometryType = 'Surface';
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = SurfaceMesh(cParams)
            obj = obj@Mesh(cParams);
            obj.init(cParams)
            
        end
        
        function hMin = computeMinCellSize(obj)
            x1 = obj.coord(obj.connec(:,1),:);
            x2 = obj.coord(obj.connec(:,2),:);
            x3 = obj.coord(obj.connec(:,3),:);
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
            x1 = obj.coord(obj.connec(:,1),:);
            x2 = obj.coord(obj.connec(:,2),:);
            x3 = obj.coord(obj.connec(:,3),:);
            x1x2 = (x2-x1);
            x2x3 = (x3-x2);
            x1x3 = (x1-x3);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            n23 = sqrt(x2x3(:,1).^2 + x2x3(:,2).^2);
            n13 = sqrt(x1x3(:,1).^2 + x1x3(:,2).^2);
            hs = max([n12,n23,n13],[],2);
            hMean = max(hs);
        end
        
        function plot(obj) %Black
            faceColor = "red";
            faceAlpha = 0.3;
            edgeAlpha = 0.5;
            if size(obj.connec,2) == 3 && size(obj.coord,2) == 3
                x = obj.coord(:,1);
                y = obj.coord(:,2);
                z = obj.coord(:,3);
                p = trisurf(obj.connec,x,y,z);
                p.FaceColor = [1 0 0];
                p.FaceAlpha = 1;
                p.EdgeColor = 'none';
                hold on
            else
                p = patch('vertices',obj.coord,'faces',obj.connec);
                p.EdgeAlpha = edgeAlpha;
                p.EdgeLighting = 'flat';
                p.FaceColor = faceColor;%[167,238,237]/265; 'green';'red';%
                p.FaceLighting = 'flat';
                p.FaceAlpha = faceAlpha;
                p.LineWidth = 1.5;
                axis('equal');
                hold on
            end
        end

        function m = provideExtrudedMesh(obj, height)
            s.unfittedMesh = obj;
            s.height       = height;
            me = MeshExtruder(s);
            m = me.extrude();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end