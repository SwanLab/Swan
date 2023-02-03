classdef FullInnerMeshCreator_Matlab < InnerMeshExporter
    
    methods (Access = public)
        
        function obj = FullInnerMeshCreator_Matlab(cParams)
            obj.init(cParams)
            
        end

        function m = export(obj)
            coordInner     = obj.unfittedMesh.innerMesh.mesh.coord;
            connecInner    = obj.unfittedMesh.innerMesh.mesh.connec;
            coordCutInner  = obj.unfittedMesh.innerCutMesh.mesh.coord;
            connecCutInner = obj.unfittedMesh.innerCutMesh.mesh.connec;
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