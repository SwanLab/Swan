classdef FullInnerMeshCreator_GiD2 < FullInnerMeshCreator
    properties (Access = private)
        filename
        meshFileName
        meshElementSize
        swanPath
        gidPath
        tclPath
        resFilePath
    end
    
    methods (Access = public)
        
        function obj = FullInnerMeshCreator_GiD2(cParams)
            obj.init(cParams)
        end

        function m = export(obj)
            obj.exportMshThroughGiD();
            m = obj.readMsh();
        end

        function exportMshThroughGiD(obj)
            obj.unfittedMesh.printNew(obj.filename);
            a = 0;
            s2g = SwanGiDInterface(a);
            resFile = obj.getResFilePath();
            s2g.generateMesh(resFile);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.unfittedMesh    = cParams.unfittedMesh;
            obj.filename        = cParams.filename;
            obj.meshElementSize = cParams.meshElementSize;
            obj.meshFileName    = cParams.meshFileName;
            obj.swanPath        = cParams.swanPath;
            obj.gidPath         = cParams.gidPath;
            obj.tclPath         = [obj.swanPath,'PostProcess/STL/'];
            obj.resFilePath     = obj.getResFilePath();
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