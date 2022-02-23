classdef MeshPlotter < handle
    
    properties (Access = private)
        mesh
        isBackground
        faceColor
        faceAlpha
        edgeAlpha
    end
      
    methods (Access = public)
        
        function obj = MeshPlotter(cParams)
            obj.init(cParams);
        end

        function plot(obj)
            m = obj.mesh;
            if obj.isBackground
                obj.plotBackgroundMesh();
            else
                if m.ndim + m.kFace == 1 
                    obj.plotLineMesh();
                elseif m.ndim + m.kFace == 2
                    obj.plotSurfaceMesh();
                end
                
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh         = cParams.mesh;
            obj.isBackground = cParams.isBackground;
            obj.faceColor    = cParams.faceColor;
            obj.faceAlpha    = cParams.faceAlpha;
            obj.edgeAlpha    = cParams.edgeAlpha;
        end
        
        function plotBackgroundMesh(obj)
            coor = obj.mesh.coord;
            conn = obj.mesh.connec;
            p = patch('vertices',coor,'faces',conn);
            p.EdgeColor = 'k';
            p.EdgeAlpha = 0.3;
            p.EdgeLighting = 'flat';
            p.FaceColor = 'none';
            p.FaceLighting = 'flat';
            p.FaceAlpha = 0;
            p.LineWidth = 0.5;
            axis('equal');
        end
        
        function plotLineMesh(obj)
            m = obj.mesh;
            p = patch('vertices',m.coord,'faces',m.connec);
            p.EdgeColor = 'b';
            p.EdgeAlpha = 1;
            p.EdgeLighting = 'flat';
            p.LineWidth = 0.5;
            p.LineStyle = '-';
            axis('equal');
            nodes = unique(m.connec(:));
            if size(m.coord,2) == 3
                x = m.coord(:,1);
                y = m.coord(:,2);
                z = m.coord(:,3);
                hold on
                p = plot3(x(nodes),y(nodes),z(nodes),'.r');
                p.MarkerSize = 6;
            else
                x = m.coord(:,1);
                y = m.coord(:,2);
                hold on
                p = plot(x(nodes),y(nodes),'.r');
                p.MarkerSize = 14;
            end
        end
        
        
%         function plotSurfaceMesh(obj)
%             m = obj.mesh;
%             if size(m.connec,2) == 3 && size(m.coord,2) == 3
%                 x = m.coord(:,1);
%                 y = m.coord(:,2);
%                 z = m.coord(:,3);
%                 p = trisurf(m.connec,x,y,z);
%                 p.FaceColor = [1 0 0];
%                 p.FaceAlpha = 1;
%                 p.EdgeColor = 'none';
%                 hold on
%             else
%                 p = patch('vertices',m.coord,'faces',m.connec);
%                 p.EdgeAlpha = 0.5;
%                 p.EdgeLighting = 'flat';
%                 p.FaceColor = 'red';%[167,238,237]/265; 'green';'red';%
%                 p.FaceLighting = 'flat';
%                 p.FaceAlpha = 0.3;
%                 p.LineWidth = 1.5;
%                 axis('equal');
%                 hold on
%             end
%             
%         end
        
        function plotSurfaceMesh(obj) %Black
            m = obj.mesh;
            if size(m.connec,2) == 3 && size(m.coord,2) == 3
                x = m.coord(:,1);
                y = m.coord(:,2);
                z = m.coord(:,3);
                p = trisurf(m.connec,x,y,z);
                p.FaceColor = [1 0 0];
                p.FaceAlpha = 1;
                p.EdgeColor = 'none';
                hold on
            else
                p = patch('vertices',m.coord,'faces',m.connec);
                p.EdgeAlpha = obj.edgeAlpha;
                p.EdgeLighting = 'flat';
                p.FaceColor = obj.faceColor;%[167,238,237]/265; 'green';'red';%
                p.FaceLighting = 'flat';
                p.FaceAlpha = obj.faceAlpha;
                p.LineWidth = 1.5;
                axis('equal');
                hold on
            end
            
        end
        
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = MeshPlotterFactory();
            obj = f.create(cParams);
        end
        
    end
    
    
end

