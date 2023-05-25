classdef MeshExtruder < handle
    properties (Access = private)
        filename
        swanPath
        gidPath
        tclPath
        resFilePath
    end
    
    methods (Access = public)
        
        function obj = MeshExtruder(cParams)
            obj.init(cParams)
        end

        function extrude(obj)
            obj.writeTCLfile();
            obj.writeExportTCLfile();
            obj.runGiD()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filename        = cParams.filename;
            obj.swanPath        = '/home/ton/Github/Swan/';
            obj.gidPath         = '/home/ton/GiDx64/gid-16.1.2d/';
            obj.tclPath         = [obj.swanPath,'PostProcess/STL/'];
%             obj.resFilePath     = obj.getResFilePath();
            obj.resFilePath     = '/home/ton/Github/Swan/Output/hellothere/hellothere1.flavia.res';
        end

        function writeTCLfile(obj)
            tclFile = [obj.tclPath,'callGiD.tcl'];
            stlFileTocall = 'CreateSurfaceSTL.tcl';
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
            fid = fopen(tclFile,'w+');
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set output "$path/sampleMesh" \n']);
            fprintf(fid,['set inputFile "',obj.resFilePath,'"\n']);
            fprintf(fid,['set meshFile "$path/sampleMesh" \n']);
            fprintf(fid,['set gidProjectName "$path/sampleMesh" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['CreateSurfaceSTL $inputFile $output $meshFile $gidProjectName $gidBasPath \n']);
            fclose(fid);
        end

        function writeExportTCLfile(obj)
            tclFile = [obj.tclPath,'.tcl'];
            stlFileTocall = 'ExportSTL.tcl';
            fid = fopen(tclFile,'w+');
            gidBasPath = [obj.gidPath,'templates/STL.bas'];
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set input "$path/sampleMesh.gid" \n']);
            fprintf(fid,['set output "$path/sampleMeshFile.stl" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['ExportSTL $input $output $gidBasPath \n']);
            fclose(fid);
        end

        function runGiD(obj)
            command = [obj.gidPath,'gid -t "source ',obj.tclPath,'callGiD.tcl"'];
            unix(command);
            
            command = [obj.gidPath,'gid -t "source ',obj.tclPath,'callGiD2.tcl"'];
            unix(command);
        end
        
    end
    
end