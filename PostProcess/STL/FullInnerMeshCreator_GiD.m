classdef FullInnerMeshCreator_GiD < FullInnerMeshCreator
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
        
        function obj = FullInnerMeshCreator_GiD(cParams)
            obj.init(cParams)
        end

        function m = export(obj)
            obj.exportMshThroughGiD();
            m = obj.readMsh();
        end

        function exportMshThroughGiD(obj)
            obj.unfittedMesh.print(obj.filename);
            obj.createSurfaceMeshTclFromTemplate();
            tclFile = [obj.tclPath,'CreateSurfaceMeshFile.tcl" '];
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',tclFile];
            system(command);
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
            f = [obj.gidPath, obj.meshFileName,'.msh'];
        end

        function f = getResFilePath(obj)
            name = obj.filename;
            swan = obj.swanPath;
            f = [swan, 'Output/', name , '/', name, '1.flavia.res'];
        end

        function createSurfaceMeshTclFromTemplate(obj)
            templateText = fileread('CreateSurfaceMeshFile_Template.tcl');
            targetFile = ['PostProcess/STL/', 'CreateSurfaceMeshFile.tcl'];
            fid = fopen(targetFile, 'w');
            fprintf(fid, ['set input_post_res "', obj.resFilePath, '"\n']);
            fprintf(fid, ['set mesh_element_size "', obj.meshElementSize, '"\n']);
            fprintf(fid, ['set mesh_name "', obj.meshFileName, '"\n']);
            fprintf(fid, ['set gidpath "', obj.gidPath, '"\n']);
            fprintf(fid, ['\n',templateText]);
            fclose(fid);
        end
        
    end
    
end