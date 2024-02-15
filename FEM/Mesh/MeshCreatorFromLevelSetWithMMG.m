classdef MeshCreatorFromLevelSetWithMMG < handle
    
    properties (Access = private)
        newMesh
        meshBackgroundMmg
        unfittedMeshMmg
    end
    
    properties (Access = private)
        meshBackground
        levelSetValue
        fileName
        hausdV
        hMinV
        hMaxV
        hMesh
    end
    
    methods (Access = public)
        
        function obj = MeshCreatorFromLevelSetWithMMG(cParams)
            obj.init(cParams)
        end

        function m = create(obj)
            obj.createMeshBackgroundMmg();
            obj.setMeshBackgroundMmgParameters();
            obj.createUnfittedMeshMmg();
            obj.createMesh();
            obj.plotMesh();
            m = obj.newMesh;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.meshBackground = cParams.meshBackground;
            obj.levelSetValue = cParams.levelSet;
            obj.fileName = cParams.fileName;
            obj.hMesh = cParams.hMesh;
            %obj.hausdV = 0.01;%0.005
            %obj.hMinV  = 0.001;%0.0005
            %obj.hMaxV  = 0.01;%0.005
            obj.hausdV = obj.hMesh;%0.005
            obj.hMinV  = obj.hMesh/10;%0.0005
            obj.hMaxV  = obj.hMesh;%0.005
        end
        
        function createMeshBackgroundMmg(obj)
            mshG = obj.createMeshBackgroundGypsiLab();
            obj.meshBackgroundMmg  = mmg(mshG,1e-3);
        end
        
        function mshG = createMeshBackgroundGypsiLab(obj)
            mesh  = obj.meshBackground;
            coord = mesh.coord;
            coord(:,3) = 0;
            mshG = msh(coord,mesh.connec);
        end
        
        function setMeshBackgroundMmgParameters(obj) 
            mesh = obj.meshBackgroundMmg;
            hausd(mesh,obj.hausdV);
            hmin(mesh,obj.hMinV);
            hmax(mesh,obj.hMaxV);
            mesh.oldFileName = [obj.fileName,'Out'];
            mesh.newFileName = [obj.fileName,'In'];
            map(mesh,obj.levelSetValue);
        end
        
        function createUnfittedMeshMmg(obj)
            uMesh = runLs(obj.meshBackgroundMmg);
            obj.unfittedMeshMmg = uMesh;
        end
        
        function m = createMesh(obj)
            [s.coord,s.connec] = obj.obtainCoordsAndConnec();
            obj.newMesh = Mesh.create(s);
        end
        
        function [coord,connec] = obtainCoordsAndConnec(obj)
            uMesh = obj.unfittedMeshMmg;
            it = uMesh.col == 3;
            connec = uMesh.elt(it,:);
            coord  = uMesh.vtx(:,1:2);
            [coord,connec] = obj.computeUniqueCoordConnec(coord,connec);
        end
        
        function plotMesh(obj)
            figure(100)
            clf
            obj.newMesh.plot();
            drawnow
        end
        
    end
    
    methods (Access = private, Static)
        
        function [newCoord,newConnec] = computeUniqueCoordConnec(coord,connec)
            allNodes = connec(:);
            [uNodes,ind,ind2] = unique(allNodes,'rows','stable');
            allCoords    = coord;
            uniqueCoords = allCoords(uNodes,:);
            newCoord    = uniqueCoords;
            
            nnode = size(connec,2);
            nCell = size(connec,1);
            newConnec = reshape(ind2,nCell,nnode);
        end

    end
    
end