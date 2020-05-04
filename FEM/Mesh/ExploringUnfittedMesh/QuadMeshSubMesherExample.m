classdef QuadMeshSubMesherExample < MeshAndGaussOfUnfittedQuadMesh
    
    methods (Access = public)
        
        function obj = QuadMeshSubMesherExample()
            
            obj.init();
            close all
%             coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
%             connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];
            nameCase = 'cutMeshProvisionalQuad';
            
            [connecFull,connecCut,cutElems] = obj.computeConnecCutAndFull(obj.levelSet,obj.backgroundMesh.connec);
            mInterior = obj.computeMesh(connecFull,obj.backgroundMesh.coord);
            
            sM.coord  = obj.backgroundMesh.coord;
            sM.connec = connecCut;
            backgroundCutMesh = Mesh().create(sM);
            
            s.backgroundMesh = backgroundCutMesh;
            s.cutElems = cutElems;
            s.levelSet  = obj.levelSet;
            s.lastNode = max(obj.backgroundMesh.connec(:));
            cutMesh = CutMeshProvisionalQuadrilater(s);
            cutMesh.compute();
            
            mCutInterior = obj.computeMesh(cutMesh.connec,cutMesh.coord);
            figure
            mInterior.plot();
            hold on
            mCutInterior.plot();
            
            d = load(nameCase);
            
            n1 = norm(d.connec(:) - cutMesh.connec(:));
            n2 = norm(d.xCoordsIso(:) - cutMesh.xCoordsIso(:));
            n3 = norm(d.cellContainingSubcell - cutMesh.cellContainingSubcell);
            n4 = norm(d.coord(:) - cutMesh.coord(:));
            
            error = n1 + n2 + n3 + n4
            
        end
        
    end
    
    
end