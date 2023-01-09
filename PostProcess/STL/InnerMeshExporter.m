classdef InnerMeshExporter < handle
    
    % Rename to InnerMeshExporter(FromUnfittedMesh)
    % Two types:
    %   - matlab: 
    %   - gid: we will NOT use it for exporting to STL, just for getting
    %   the inner mesh.

    properties (Access = public)
        
    end
    
    properties (Access = private) % Inputs
        filename
        type
        unfittedMesh
    end
    
    properties (Access = private) % Calculated
    end
    
    methods (Access = public)
        
        function obj = InnerMeshExporter(cParams)
            obj.init(cParams)
        end

        function m = export(obj)
            switch obj.type
                case 'Matlab'
                    m = obj.exportUsingMatlab();
                case 'GiD'
                    m = obj.exportUsingGiD();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.unfittedMesh  = cParams.unfittedMesh;
            obj.type          = cParams.type;
            obj.filename      = cParams.filename;% only for gid
        end

        function m = exportUsingMatlab(obj)
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

        function m = exportUsingGiD(obj)
            mshFile = obj.unfittedMesh.exportGiD();
            m = obj.readMsh(mshFile);
            % return msh in Mesh()
        end

        function m = readMsh(obj, file)
        end
        
    end
    
end