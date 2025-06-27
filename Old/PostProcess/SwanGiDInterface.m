classdef SwanGiDInterface < handle
    
    properties (Access = public)
        swanPath
        gidPath
        tclPath

        gidMode = 'gid_offscreen -offscreen '; % 'gid '
%         gidMode = 'gid '; % 'gid '
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = SwanGiDInterface()
            run('UserVariables.m')
            obj.swanPath = pwd;
            obj.gidPath  = gid_path;
            obj.tclPath  = fullfile([obj.swanPath, '/PostProcess/STL/']);
        end

        function generateMesh(obj, resultsFile)
            obj.writeSurfaceTclFile(resultsFile);
            obj.writeGenerateMeshTclFile();
            obj.runSurfaceTcl();
            obj.runGenerateMeshTcl();
            obj.cleanupGenerateMesh();
        end

        function extrudeMesh(obj, mesh, height)
            obj.writeExportMshTclFile(mesh);
            obj.runExportMshTcl();
            obj.writeExtrudeTclFile(height);
            obj.writeGenerateMeshTclFile();
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

    end

    methods (Access = private)

        function writeSurfaceTclFile(obj, resultsFile)
            tclFile = fullfile([obj.tclPath,'callGiD_CreateSurface.tcl']);
            stlFileTocall = 'CreateSurfaceNew.tcl';
            gidBasPath = fullfile([obj.gidPath,'templates\DXF.bas']);
            fid = fopen(tclFile,'w+');
            fprintf(fid,['set path "',obj.tclPath,'"\n']);
            fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
            fprintf(fid,['source $path$tclFile \n']);
            fprintf(fid,['set output "$path\sampleMesh" \n']);
            fprintf(fid,['set inputFile "',resultsFile,'"\n']);
            fprintf(fid,['set meshFile "$path\sampleMesh" \n']);
            fprintf(fid,['set gidProjectName "$path\sampleMesh" \n']);
            fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
            fprintf(fid,['CreateSurfaceNew $inputFile $output $meshFile $gidProjectName $gidBasPath \n']);

%             fprintf(fid, 'set path "%s"\n', obj.tclPath);
%             fprintf(fid, 'set tclFile "%s"\n', stlFileTocall);
%             fprintf(fid, 'source $path$tclFile\n');
%             fprintf(fid, 'set output "$path\sampleMesh"\n');
%             fprintf(fid, 'set inputFile "%s"\n', resultsFile);
%             fprintf(fid, 'set meshFile "$path\sampleMesh"\n');
%             fprintf(fid, 'set gidProjectName "$path\sampleMesh"\n');
%             fprintf(fid, 'set gidBasPath "%s"\n', gidBasPath);
%             fprintf(fid, 'CreateSurfaceNew $inputFile $output $meshFile $gidProjectName $gidBasPath\n');

            fclose(fid);
        end

        function runSurfaceTcl(obj)
            % Create Surface
%             command = obj.callGiDFunctionFile('callGiD_CreateSurface.tcl"');
            command = ['"', ...
            'C:\Program Files\GiD\GiD 17.0.4\gid_offscreen', ...
            '" -offscreen -t "source ', ...
            '"C:/Users/Biel/Desktop/UNI/TFG/PostProcess/STL/callGiD_CreateSurface.tcl"', ...
            '"'];
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
            fprintf(fid,['set gidProjectName "$path\sampleMesh" \n']);
            fprintf(fid,['GenerateMesh $gidProjectName \n']);

%             fprintf(fid, 'set path "%s"\n', obj.tclPath);
%             fprintf(fid, 'set tclFile "%s"\n', stlFileTocall);
%             fprintf(fid, 'source $path$tclFile\n');
%             fprintf(fid, 'set gidProjectName "$path\sampleMesh"\n');
%             fprintf(fid, 'GenerateMesh $gidProjectName\n');
            fclose(fid);
        end

        function runGenerateMeshTcl(obj)
            % Generate Mesh
            % command = obj.callGiDFunctionFile('callGiD_GenerateMesh.tcl"');
            command = ['"', ...
                'C:\Program Files\GiD\GiD 17.0.4\gid_offscreen', ...
                '" -offscreen -t "source ', ...
                '"C:/Users/Biel/Desktop/UNI/TFG/PostProcess/STL/callGiD_GenerateMesh.tcl"', ...
                '"'];
            system(command);
        end

        function writeExtrudeTclFile(obj, height)
            tclFile = [obj.tclPath,'callGiD_Extrude.tcl'];
            stlFileTocall = 'ExtrudeSurface.tcl';
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
            fid = fopen(tclFile,'w+');
%             fprintf(fid,['set path "',obj.tclPath,'"\n']);
%             fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
%             fprintf(fid,['set mshname "sampleMesh.msh" \n']);
%             fprintf(fid,['source $path$tclFile \n']);
%             fprintf(fid,['set input "$path/HmmLetMeCook.msh" \n']);
%             fprintf(fid,['set output "$path/sampleMesh" \n']);
%             fprintf(fid,['set height "', sprintf('%f', height),'" \n']);
%             fprintf(fid,['ExtrudeSurface $input $height $output \n']);
            fprintf(fid, 'set path "%s"\n', obj.tclPath);
            fprintf(fid, 'set tclFile "%s"\n', stlFileTocall);
            fprintf(fid, 'set mshname "sampleMesh.msh"\n');
            fprintf(fid, 'source $path$tclFile\n');
            fprintf(fid, 'set input "$path/HmmLetMeCook.msh"\n');
            fprintf(fid, 'set output "$path/sampleMesh"\n');
            fprintf(fid, 'set height "%f"\n', height);
            fprintf(fid, 'ExtrudeSurface $input $height $output\n');
            fclose(fid);
        end

        function runExtrudeTcl(obj)
            % Extrude Surface
            command = obj.callGiDFunctionFile('callGiD_Extrude.tcl"');
            system(command);
        end

        function writeExportSTLTclFile(obj)
            tclFile = [obj.tclPath,'callGiD_ExportSTL.tcl'];
            stlFileTocall = 'ExportSTL.tcl';
            fid = fopen(tclFile,'w+');
            gidBasPath = [obj.gidPath,'templates/STL.bas'];
%             fprintf(fid,['set path "',obj.tclPath,'"\n']);
%             fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
%             fprintf(fid,['source $path$tclFile \n']);
%             fprintf(fid,['set input "$path/HmmLetMeCook.msh" \n']);
%             fprintf(fid,['set output "$path/sampleMeshFile.stl" \n']);
%             fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
%             fprintf(fid,['ExportSTL $input $output $gidBasPath \n']);
            fprintf(fid, 'set path "%s"\n', obj.tclPath);
            fprintf(fid, 'set tclFile "%s"\n', stlFileTocall);
            fprintf(fid, 'source $path$tclFile\n');
            fprintf(fid, 'set input "$path/HmmLetMeCook.msh"\n');
            fprintf(fid, 'set output "$path/sampleMeshFile.stl"\n');
            fprintf(fid, 'set gidBasPath "%s"\n', gidBasPath);
            fprintf(fid, 'ExportSTL $input $output $gidBasPath\n');
            fclose(fid);
        end

        function runExportSTLTcl(obj)
            % Export STL
            command = obj.callGiDFunctionFile('callGiD_ExportSTL.tcl"');
            system(command);
        end

        function writeExportMshTclFile(obj, mesh)
            % Print mesh
            mesh.print('TempMeshFile', 'GiD');

            tclFile = [obj.tclPath,'callGiD_ExportMSH.tcl'];
            stlFileTocall = 'ExportMSH.tcl';
            fid = fopen(tclFile,'w+');
            gidBasPath = [obj.gidPath,'templates/DXF.bas'];
%             fprintf(fid,['set path "',obj.tclPath,'"\n']);
%             fprintf(fid,['set swanpath "',obj.swanPath,'"\n']);
%             fprintf(fid,['set tclFile "',stlFileTocall,'"\n']);
%             fprintf(fid,['source $path$tclFile \n']);
%             fprintf(fid,['set input "$swanpath/TempMeshFile.flavia.res" \n']);
% %             fprintf(fid,['set output "$swanpath/resultingMesh.stl" \n']);
%             fprintf(fid,['set gidBasPath "',gidBasPath,'" \n']);
%             fprintf(fid,['ExportMSH $input $gidBasPath \n']);
            fprintf(fid, 'set path "%s"\n', obj.tclPath);
            fprintf(fid, 'set swanpath "%s"\n', obj.swanPath);
            fprintf(fid, 'set tclFile "%s"\n', stlFileTocall);
            fprintf(fid, 'source $path$tclFile\n');
            fprintf(fid, 'set input "$swanpath/TempMeshFile.flavia.res"\n');
            fprintf(fid, 'set gidBasPath "%s"\n', gidBasPath);
            fprintf(fid, 'ExportMSH $input $gidBasPath\n');
            fclose(fid);
        end

        function runExportMshTcl(obj)
            % Export MSH
            command = obj.callGiDFunctionFile('callGiD_ExportMSH.tcl"');
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
            obj.cleanupExportSTL();
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

        function c = callGiDFunctionFile(obj, file)
            c = [obj.gidPath,obj.gidMode,'-t "source ',obj.tclPath, file];
        end

    end

end

% command = ['"', ...
%     'C:\Program Files\GiD\GiD 17.0.4\gid_offscreen', ...
%     '" -offscreen -t "source ', ...
%     '"C:/Users/Biel/Desktop/UNI/TFG/PostProcess/STL/callGiD_GenerateMesh.tcl"', ...
%     '"'];

% command = ['"', ...
%     'C:\Program Files\GiD\GiD 17.0.4\gid_offscreen', ...
%     '" -offscreen -t "source ', ...
%     '"C:/Users/Biel/Desktop/UNI/TFG/PostProcess/STL/callGiD_CreateSurface.tcl"', ...
%     '"'];
