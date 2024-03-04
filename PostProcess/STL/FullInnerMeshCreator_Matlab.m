classdef FullInnerMeshCreator_Matlab < FullInnerMeshCreator
    
    methods (Access = public)
        
        function obj = FullInnerMeshCreator_Matlab(cParams)
            obj.init(cParams)
            
        end

        function m = export(obj)
            uM = obj.unfittedMesh;
            innerMesh    = uM.innerMesh.mesh;
            innerCutMesh = uM.innerCutMesh.mesh;
            coordCutInner  = innerCutMesh.coord;
            connecCutInner = innerCutMesh.connec;

            switch uM.innerMesh.mesh.type
                case 'TRIANGLE'
                    coordInner     = innerMesh.coord;
                    connecInner    = innerMesh.connec;

                case 'QUAD'
                    innerMeshQuad = innerMesh;
                    q2t = QuadToTriMeshConverter();
                    innerMeshTri = q2t.convert(innerMeshQuad, innerMeshQuad.nnodes);
                    innerMeshTri = innerMeshTri.computeCanonicalMesh();
                    coordInner  = innerMeshTri.coord;
                    connecInner = innerMeshTri.connec;
            end
            ncoord = size(coordInner,1);
            % Find duplicate coordinates
            [~,oldNodes,newNodes] = intersect(coordCutInner,coordInner,'rows');
            old2new = (1:1:max(max(connecCutInner))) + ncoord;
            old2new(oldNodes) = newNodes;
            newConnecCutInner = old2new(connecCutInner);
            s.coord  = [coordInner; coordCutInner];
            s.connec = [connecInner; newConnecCutInner];
            msh = Mesh.create(s);
            m = msh.computeCanonicalMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.unfittedMesh  = cParams.unfittedMesh;
        end
        
    end
    
end