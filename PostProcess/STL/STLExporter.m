classdef STLExporter < handle
    properties (Access = private)
        filename
        mesh
    end
    
    methods (Access = public)
        
        function obj = STLExporter(cParams)
            obj.init(cParams)
        end

        function m = export(obj)
            a = 0;
            s2g = SwanGiDInterface(a);
            s2g.exportSTL(obj.mesh);
            m = obj.readMsh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filename = cParams.filename;
            obj.mesh     = cParams.mesh;
        end

        function m = readMsh(obj)
            s.filePath = obj.getOutputFileName();
            mR = MshReader(s);
            m = mR.read();
        end

        function f = getOutputFileName(obj)
%             f = [obj.gidPath, obj.meshFileName,'.msh'];
            f = '/home/ton/Github/Swan/PostProcess/STL/sampleMesh.msh';
        end

        function f = getResFilePath(obj)
            name = obj.filename;
            swan = obj.swanPath;
            f = [swan, name, '.flavia.res'];
        end
        
    end
    
end