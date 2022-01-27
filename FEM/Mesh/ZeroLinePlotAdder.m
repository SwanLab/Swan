classdef ZeroLinePlotAdder < handle
    
    properties (Access = private)
        x
        y
        z
        figID
        meshBackground
        unfittedMesh
    end
    
    methods (Access = public)
        
        function obj = ZeroLinePlotAdder(cParams)
            obj.init(cParams);
        end
        
        function addLine(obj)
            obj.createMeshBackground();
            obj.createUnfittedMesh();
            obj.plotLine();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.x = cParams.x;
            obj.y = cParams.y;
            obj.z = cParams.z;
            obj.figID = cParams.figID;
        end
        
        function createMeshBackground(obj)
            [xp,yp] = meshgrid(obj.x,obj.y);
            xp = xp(:);
            yp = yp(:);
            tri = delaunay(xp,yp);
            s.coord = [xp yp];
            s.connec = tri;
            mesh = Mesh_Total(s);
            obj.meshBackground = mesh;
        end
        
        function createUnfittedMesh(obj)
            s.meshBackground = obj.meshBackground;
            s.backgroundMesh = obj.meshBackground.innerMeshOLD;
            s.boundaryMesh   = obj.meshBackground.boxFaceMeshes;
            s.unfittedType   = 'BOUNDARY';
            s.isInBoundary = false;
            cParams = SettingsMeshUnfitted(s);
            uMesh = UnfittedMesh(cParams);
            levelSet = obj.z;
            uMesh.compute(levelSet);
            obj.unfittedMesh = uMesh;
        end
        
        function plotLine(obj)
            figure(obj.figID)
            patch('vertices',obj.unfittedMesh.coord,...
                  'faces',obj.unfittedMesh.connec,...
                  'edgecolor',[0.5 0 0],...
                   'edgealpha',0.5,...
                   'edgelighting','flat',...
                   'facecolor',[1 0 0],...
                   'facelighting','flat',...
                   'LineWidth',6)
            axis('equal')
        end
        
    end

end