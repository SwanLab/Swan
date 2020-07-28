classdef MeshPlotter < handle
    
    properties (Access = private)
        mesh
        isBackground
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
        end
        
        function plotBackgroundMesh(obj)
            coor = obj.mesh.coord;
            conn = obj.mesh.connec;
            p = patch('vertices',coor,'faces',conn);
            p.EdgeColor = 'k';
            p.EdgeAlpha = 0.1;
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
            p.LineWidth = 2.5;
            p.LineStyle = '-';
            axis('equal');
            nodes = unique(m.connec(:));            
            if size(m.coord,2) == 3
                x = m.coord(:,1);
                y = m.coord(:,2);
                z = m.coord(:,3);                   
                hold on
                plot3(x(nodes),y(nodes),z(nodes),'.r')                
            else
                x = m.coord(:,1);
                y = m.coord(:,2);
                hold on
                plot(x(nodes),y(nodes),'.r')                                
            end            
        end
        
        function plotSurfaceMesh(obj)
            m = obj.mesh;
            nodes = unique(m.connec(:));            
            if size(m.connec,2) == 3 && size(m.coord,2) == 3
                x = m.coord(:,1);
                y = m.coord(:,2);
                z = m.coord(:,3);
                p = trisurf(m.connec,x,y,z);
                p.FaceColor = 'cyan';
                p.FaceAlpha = 0.8;
                axis equal;
                hold on
                plot3(x(nodes),y(nodes),z(nodes),'.r')
            else
                x = m.coord(:,1);
                y = m.coord(:,2);                
                p = patch('vertices',m.coord,'faces',m.connec);
                p.EdgeAlpha = 0.5;
                p.EdgeLighting = 'flat';
%                p.FaceColor = [1 0 0];
                p.FaceColor = 'cyan';                
                p.FaceLighting = 'flat';
                p.FaceAlpha = 1;
                p.LineWidth = 1.5;
                axis('equal');
                hold on
                plot(x(nodes),y(nodes),'.r')                
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

