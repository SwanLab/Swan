classdef SwanGiDInterface < handle
    
    properties (Access = public)
        swanPath
        gidPath
        tclPath
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = SwanGiDInterface(cParams)
            obj.swanPath = pwd;
            obj.gidPath  = '/home/ton/GiDx64/gid-16.1.2d/';
            obj.tclPath  = [obj.swanPath, '/PostProcess/STL/'];
        end

        function generateMesh(obj, resultsFile)
            obj.writeSurfaceTclFile(resultsFile);
            obj.writeGenerateMeshTclFile();
            obj.runSurfaceTcl();
            obj.runGenerateMeshTcl();
            obj.cleanupGenerateMesh();
        end

        function extrudeMesh(obj, mesh)
            obj.writeExportMshTclFile(mesh);
            obj.runExportMshTcl();
%             obj.writeSurfaceTclFile(resultsFile);
            obj.writeExtrudeTclFile();
            obj.writeGenerateMeshTclFile();
%             obj.runSurfaceTcl();
            obj.runExtrudeTcl();
            obj.runGenerateMeshTcl();
            obj.cleanupExtrudeMesh();
        end
        
        function exportSTL(obj, mesh)
            obj.writeExportMshTclFile(mesh);
            obj.runExportMshTcl();
            obj.writeExportSTLTclFile();
            obj.runExportSTLTcl();
            obj.cleanupExportSTL();
        end

        function extrudeAndExport(obj)
            resultsFile = '/home/ton/Github/Swan/hellouNou.flavia.res';
            obj.writeSurfaceTclFile(resultsFile);
            obj.writeExtrudeTclFile();
            obj.writeGenerateMeshTclFile();
            obj.writeExportSTLTclFile();
            obj.runSurfaceTcl();
            obj.runExtrudeTcl();
            obj.runGenerateMeshTcl();
            obj.runExportSTLTcl();
        end

    end

    methods (Access = private)

        function writeSurfaceTclFile(obj, resultsFile)
            tclFile = [obj.tclPath,'callGiD_CreateSurface.tcl'];
            stlFileTocall = 'CreateSurfaceNew.tcl';
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
            fid = fopen(tclFile,'w+');
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set output "$path/sampleMesh" \n']);
            fprintf(fid,['set inputFile "',resultsFile,'"\n']);
            fprintf(fid,['set meshFile "$path/sampleMesh" \n']);
            fprintf(fid,['set gidProjectName "$path/sampleMesh" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['CreateSurfaceNew $inputFile $output $meshFile $gidProjectName $gidBasPath \n']);
            fclose(fid);
        end

        function runSurfaceTcl(obj)
            % Create Surface
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',obj.tclPath,'callGiD_CreateSurface.tcl"'];
            system(command);
        end

        function writeGenerateMeshTclFile(obj)
            tclFile = [obj.tclPath,'callGiD_GenerateMesh.tcl'];
            stlFileTocall = 'GenerateMesh.tcl';
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
            fid = fopen(tclFile,'w+');
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set gidProjectName "$path/sampleMesh" \n']);
            fprintf(fid,['GenerateMesh $gidProjectName \n']);
            fclose(fid);
        end

        function runGenerateMeshTcl(obj)
            % Generate Mesh
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',obj.tclPath,'callGiD_GenerateMesh.tcl"'];
            system(command);
        end

        function writeExtrudeTclFile(obj)
            tclFile = [obj.tclPath,'callGiD_Extrude.tcl'];
            stlFileTocall = 'ExtrudeSurface.tcl';
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
            fid = fopen(tclFile,'w+');
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['set mshname "sampleMesh.msh" \n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set input "$path/HmmLetMeCook.msh" \n']);
            fprintf(fid,['set output "$path/sampleMesh" \n']);
            fprintf(fid,['set height "0.16" \n']);
            fprintf(fid,['ExtrudeSurface $input $height $output \n']);
            fclose(fid);
        end

        function runExtrudeTcl(obj)
            % Extrude Surface
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',obj.tclPath,'callGiD_Extrude.tcl"'];
            system(command);
        end

        function writeExportSTLTclFile(obj)
            tclFile = [obj.tclPath,'callGiD_ExportSTL.tcl'];
            stlFileTocall = 'ExportSTL.tcl';
            fid = fopen(tclFile,'w+');
            gidBasPath = [obj.gidPath,'templates/STL.bas'];
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set input "$path/HmmLetMeCook.msh" \n']);
            fprintf(fid,['set output "$path/sampleMeshFile.stl" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['ExportSTL $input $output $gidBasPath \n']);
            fclose(fid);
        end

        function runExportSTLTcl(obj)
            % Export STL
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',obj.tclPath,'callGiD_ExportSTL.tcl"'];
            system(command);
        end

        function writeExportMshTclFile(obj, mesh)
            % Print mesh
            mesh.print('TempMeshFile', 'GiD');

            tclFile = [obj.tclPath,'callGiD_ExportMSH.tcl'];
            stlFileTocall = 'ExportMSH.tcl';
            fid = fopen(tclFile,'w+');
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set swanpath "',obj.swanPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set input "$swanpath/TempMeshFile.flavia.res" \n']);
%             fprintf(fid,['set output "$swanpath/resultingMesh.stl" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['ExportMSH $input $gidBasPath \n']);
            fclose(fid);
        end

        function runExportMshTcl(obj)
            % Export MSH
            command = [obj.gidPath,'gid_offscreen -offscreen -t "source ',obj.tclPath,'callGiD_ExportMSH.tcl"'];
            system(command);
            delete PostProcess/STL/HmmLetMeCook.png;
            delete PostProcess/STL/HmmLetMeCook.res;
            delete PostProcess/STL/HmmLetMeCook.vv;
        end

        % Cleanup
        function cleanupGenerateMesh(obj)
            delete PostProcess/STL/callGiD_CreateSurface.tcl
            delete PostProcess/STL/callGiD_GenerateMesh.tcl
            delete PostProcess/STL/sampleMesh
            delete PostProcess/STL/sampleMesh.png
            delete PostProcess/STL/sampleMesh.res
            delete PostProcess/STL/sampleMesh.vv
            rmdir('PostProcess/STL/sampleMesh.gid/', 's')
        end

        function cleanupExtrudeMesh(obj)
            delete PostProcess/STL/callGiD_CreateSurface.tcl
            delete PostProcess/STL/callGiD_Extrude.tcl
            delete PostProcess/STL/callGiD_GenerateMesh.tcl
            rmdir('PostProcess/STL/sampleMesh.gid/', 's')
            delete PostProcess/STL/sampleMesh
            delete PostProcess/STL/sampleMesh.png
            delete PostProcess/STL/sampleMesh.res
            delete PostProcess/STL/sampleMesh.vv
        end

        function cleanupExportSTL(obj)
            delete PostProcess/STL/callGiD_ExportMSH.tcl
            delete PostProcess/STL/callGiD_ExportSTL.tcl
            delete PostProcess/STL/HmmLetMeCook.msh
            delete TempMeshFile.flavia.msh
            delete TempMeshFile.flavia.res
        end

    end

end

