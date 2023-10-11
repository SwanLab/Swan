classdef FullInnerMeshCreator_Matlab < FullInnerMeshCreator
    
    methods (Access = public)
        
        function obj = FullInnerMeshCreator_Matlab(cParams)
            obj.init(cParams)
            
        end

        function m = export(obj)
            uM = obj.unfittedMesh;
            switch uM.innerMesh.mesh.type
                case 'TRI'
                    coordInner     = uM.innerMesh.mesh.coord;
                    connecInner    = uM.innerMesh.mesh.connec;
                    coordCutInner  = uM.innerCutMesh.mesh.coord;
                    connecCutInner = uM.innerCutMesh.mesh.connec;

                case 'QUAD'
                    innerMeshQuad = uM.innerMesh.mesh;
                    innerMeshTri = innerMeshQuad.convertToTriangleMesh();
                    coordInner  = innerMeshTri.coord;
                    connecInner = innerMeshTri.connec;
                    coordCutInner  = uM.innerCutMesh.mesh.coord;
                    connecCutInner = uM.innerCutMesh.mesh.connec;
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