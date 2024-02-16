classdef LineMesh < Mesh
    
    properties (Access = public)
        geometryType = 'Line';
        
%         coord, connec
%         kFace
    end
    
    properties (Access = private)
%         type
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LineMesh(cParams)
            obj = obj@Mesh(cParams);
            % obj.init(cParams)
        end
        
        function hMin = computeMinCellSize(obj)
            x1 = obj.coord(obj.connec(:,1),:);
            x2 = obj.coord(obj.connec(:,2),:);
            x1x2 = (x2-x1);
            n12 = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            hMin = min(n12);
        end

        function hMean = computeMeanCellSize(obj)
            x1 = obj.coord(obj.connec(:,1),:);
            x2 = obj.coord(obj.connec(:,2),:);
            x1x2 = (x2-x1);
            hs = sqrt(x1x2(:,1).^2 + x1x2(:,2).^2);
            hMean = max(hs);
        end
        
        function plot(obj)
            p = patch('vertices',obj.coord,'faces',obj.connec);
            p.EdgeColor = 'b';
            p.EdgeAlpha = 1;
            p.EdgeLighting = 'flat';
            p.LineWidth = 0.5;
            p.LineStyle = '-';
            axis('equal');
            nodes = unique(obj.connec(:));
            if size(obj.coord,2) == 3
                x = obj.coord(:,1);
                y = obj.coord(:,2);
                z = obj.coord(:,3);
                hold on
                p = plot3(x(nodes),y(nodes),z(nodes),'.r');
                p.MarkerSize = 6;
            else
                x = obj.coord(:,1);
                y = obj.coord(:,2);
                hold on
                p = plot(x(nodes),y(nodes),'.r');
                p.MarkerSize = 14;
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
        end
        
    end
    
end