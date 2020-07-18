classdef MeshPlotter < handle
    
    properties (Access = private)
        mesh
    end
    
    methods (Access = public)
        
        function obj = MeshPlotter(cParams)
            obj.init(cParams);
        end
        
        function plotWithGrid(obj)
            figure
            obj.plot();
        end
        
        function plot(obj)
            m = obj.mesh;
            if isequal(class(m),'Mesh_Total') 
                obj.plotBackgroundMesh();
            elseif isequal(class(m),'InnerMesh') 
                hold on
                if m.ndim == 2  %(InnerMesh, InnerCutMesh)
                   obj.plotInnerMeshIn2D();
                end
            else
                if ~isempty(m)
                    hold on
                    if m.ndim == 2 && m.kFace == 0 %(InnerMesh, InnerCutMesh)
                        obj.plotInnerMeshIn2D();
                    elseif m.ndim == 2 && m.kFace == -1 %(BoundaryCutMesh2D)
                        obj.plotBoundaryMesh();
                    elseif m.ndim == 3 && m.kFace == -1  %(BoundaryCutMesh3D)
                        obj.plotBoundary3DMesh();
                    end
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        
        function plotBoundary3DMesh(obj)
            m = obj.mesh;
            p = patch('vertices',m.coord,'faces',m.connec);
            %                p.Edgecolor = 'b';
            p.EdgeAlpha = 0.5;
            p.EdgeLighting = 'flat';
            p.FaceColor = [1 0 0];
            p.FaceLighting = 'flat';
            p.FaceAlpha = 1;
            p.LineWidth = 1.5;
            axis('equal');
        end
        
        function plotBoundaryMesh(obj)
            m = obj.mesh;
            p = patch('vertices',m.coord,'faces',m.connec);
            p.EdgeColor = 'g';
            p.EdgeAlpha = 1;
            p.EdgeLighting = 'flat';
            p.LineWidth = 2.5;
            p.LineStyle = '-';
            axis('equal');
        end
        
        function plotInnerMeshIn2D(obj)
            m = obj.mesh;
            p = patch('vertices',m.coord,'faces',m.connec);
            p.EdgeColor = 'b';
            p.EdgeAlpha = 1;
            p.EdgeLighting = 'flat';
            p.FaceColor = [1 0 0];
            p.FaceLighting = 'flat';
            p.LineWidth = 0.5;
            axis('equal');
        end
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
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
        
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = MeshPlotterFactory();
            obj = f.create(cParams);
        end
        
    end
    
    
end

