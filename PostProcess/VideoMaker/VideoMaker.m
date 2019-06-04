classdef VideoMaker < handle
    
    properties (GetAccess = private, SetAccess = public)
        iterations
    end
    
    properties (Access = private)
        gidPath
        tclFileName
        filesFolder
        fileName
        tclFileWriter
        fieldName
        designVariableName
    end
    
    methods (Access = public)
        
        function makeVideo(obj)
            obj.makeDesignVariableVideo();
            obj.makeRegularizedDesignVariableVideo();
        end
        
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = VideoMakerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.fileName           = cParams.caseFileName;
            obj.designVariableName = cParams.designVarType;
            obj.createPaths();
            obj.createFolder();
        end
        
    end
    
    methods (Access = private)
        
        function createFolder(obj)
            if ~exist(obj.filesFolder,'dir')
                mkdir(obj.filesFolder);
            end
        end
        
        function createPaths(obj)
            %obj.gidPath = 'C:\Program Files\GiD\GiD 13.0.4';%
            obj.gidPath = '/opt/GiDx64/13.0.2/';
            obj.filesFolder = fullfile(pwd,'Output',obj.fileName);
        end
        
        function makeDesignVariableVideo(obj)
            obj.fieldName = obj.designVariableName;
            obj.makeFieldVideo();
        end
        
        function makeRegularizedDesignVariableVideo(obj)
            obj.fieldName = 'RegularizedDensity';
            obj.makeFieldVideo();            
        end
        
        function makeFieldVideo(obj)
            obj.createTclFileName();
            obj.createTclFileWriter();
            obj.writeTclFile();
            obj.executeTclFile();
            obj.deleteTclFile();            
        end
        
        function createTclFileName(obj)
            fName = 'tcl_gid.tcl';
            obj.tclFileName = fullfile(obj.filesFolder,fName);
        end
        
        function createTclFileWriter(obj)
            cParams.type         = obj.fieldName;
            cParams.tclFileName  = obj.tclFileName;
            cParams.filesFolder  = obj.filesFolder;
            cParams.outputName   = [obj.fieldName,obj.fileName];
            cParams.iterations   = obj.iterations;
            cParams.fileName     = obj.fileName;
            obj.tclFileWriter    = TclFileWriter.create(cParams);
        end
        
        function writeTclFile(obj)
            obj.tclFileWriter.write();
        end
        
        function executeTclFile(obj)
            tFile = replace(obj.tclFileName,'\','\\');
            gFile = fullfile(obj.gidPath,'gid_offscreen');
            executingLine = ['"',gFile,'"', ' -t ' ,'"source ',tFile,'"'];
            system(executingLine);
        end
        
        function deleteTclFile(obj)
            tFile = obj.tclFileName;
            if ispc
                system(['DEL ',tFile]);
            elseif isunix
                system(['rm ',tFile]);
            elseif ismac
            end
        end
        
    end    
    
end
