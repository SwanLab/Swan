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
                    innerMeshTri = innerMeshQuad.convertToTriangleMesh();
                    innerMeshTri = innerMeshTri.computeCanonicalMesh();
                    coordInner  = innerMeshTri.coord;
                    connecInner = innerMeshTri.connec;
            end
            ncoord = size(coordInner,1);
            connecCutInner = connecCutInner + ncoord;
            s.coord  = [coordInner;  coordCutInner];
            s.connec = [connecInner; connecCutInner];
            m = Mesh(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.unfittedMesh  = cParams.unfittedMesh;
        end
        
    end
    
end