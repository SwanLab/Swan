classdef FullInnerMeshCreator_GiD < FullInnerMeshCreator
    properties (Access = private)
        filename
        gidProjectPath
        meshElementSize
        meshFileName
        swanPath
        gidPath
    end
    
    methods (Access = public)
        
        function obj = FullInnerMeshCreator_GiD(cParams)
            obj.init(cParams)
        end

        function m = export(obj)
            h = obj.unfittedMesh.innerMesh.mesh.computeMeanCellSize();
            s.filename        = 'hellothere';
            s.gidProjectPath  = '/home/ton/test_micro_project.gid';
            s.meshElementSize = num2str(h);
            s.meshFileName    = 'hmmmm22';
            s.swanPath        = '/home/ton/Github/Swan/';
            s.gidPath         = '/home/ton/GiDx64/gid-16.1.2d/';
            obj.exportMshThroughGiD(s);
            f = obj.getOutputFileName(s);
            m = obj.readMsh(f);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filename        = cParams.filename;
            obj.unfittedMesh    = cParams.unfittedMesh;
            obj.gidProjectPath  = cParams.gidProjectPath;
            obj.meshElementSize = cParams.meshElementSize;
            obj.meshFileName    = cParams.meshFileName;
            obj.swanPath        = cParams.swanPath;
            obj.gidPath         = cParams.gidPath;
        end
        
    end

    methods (Static, Access = private)

        function f = getOutputFileName(s)
            f = [s.gidPath,s.meshFileName,'.msh'];
        end

        function m = readMsh(filename)
            s.filePath = filename;
            mR = MshReader(s);
            m = mR.read();
        end
        
        function exportMshThroughGiD(obj, cParams)
            filename = cParams.filename;
            obj.print(filename);
            cParams.resFilePath = [cParams.swanPath, 'Output/',filename,'/',filename,'1.flavia.res'];
            obj.createSurfaceMeshTclFromTemplate(cParams);
            pathTcl  = [cParams.swanPath,'PostProcess/STL/'];
            
            tclFile = [pathTcl,'CreateSurfaceMeshFile.tcl" '];
            command = [cParams.gidPath,'gid_offscreen -offscreen -t "source ',tclFile];
            system(command);
        end

    end
    
end