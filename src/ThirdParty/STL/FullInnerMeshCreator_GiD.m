classdef FullInnerMeshCreator_GiD < FullInnerMeshCreator
    properties (Access = private)
        filename
        resFilePath
    end
    
    methods (Access = public)
        
        function obj = FullInnerMeshCreator_GiD(cParams)
            obj.init(cParams)
        end

        function m = export(obj)
            obj.exportMshThroughGiD();
            m = obj.readMsh();
            delete InnerMeshCreator_File.flavia.msh
            delete InnerMeshCreator_File.flavia.res
            delete src/ThirdParty/STL/sampleMesh.msh
        end

        function exportMshThroughGiD(obj)
            obj.unfittedMesh.printNew(obj.filename);
            s2g = SwanGiDInterface();
            resFile = obj.getResFilePath();
            s2g.generateMesh(resFile);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.unfittedMesh = cParams.unfittedMesh;
            obj.filename = 'InnerMeshCreator_File';
        end

        function m = readMsh(obj)
            s.filePath = obj.getOutputFileName();
            mR = MshReader(s);
            m = mR.read();
        end

        function f = getOutputFileName(obj)
            f = [pwd,'/src/ThirdParty/STL/sampleMesh.msh'];
        end

        function f = getResFilePath(obj)
            name = obj.filename;
            swan = pwd;
            f = [swan, '/', name, '.flavia.res'];
        end
        
    end
    
end