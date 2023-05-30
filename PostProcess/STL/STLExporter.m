classdef STLExporter < handle
    properties (Access = private)
        filename
        swanPath
        gidPath
        tclPath
        resFilePath
        unfittedMesh
        meshFileName
        mesh
    end
    
    methods (Access = public)
        
        function obj = STLExporter(cParams)
            obj.init(cParams)
        end

        function m =  export(obj)
%             obj.createMsh();
            obj.exportSTLThroughGiD();
            m = obj.readMsh();
        end

        function exportSTLThroughGiD(obj)
%             obj.unfittedMesh.printNew(obj.filename);
            resFile = '/home/ton/Github/Swan/PostProcess/STL/sampleMesh.msh';
            a = 0;
            s2g = SwanGiDInterface(a);
%             resFile = obj.getResFilePath();
            s2g.exportSTL(obj.mesh);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
%             obj.unfittedMesh    = cParams.unfittedMesh;
            obj.filename        = cParams.filename;
%             obj.meshElementSize = cParams.meshElementSize;
%             obj.meshFileName    = cParams.meshFileName;
%             obj.swanPath        = cParams.swanPath;
%             obj.gidPath         = cParams.gidPath;
%             obj.tclPath         = [obj.swanPath,'PostProcess/STL/'];
            obj.resFilePath     = obj.getResFilePath();
            obj.mesh = cParams.mesh;
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