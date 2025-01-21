classdef VideoMaker < handle
    
    properties (GetAccess = private, SetAccess = public)
        iterations
    end
    
    properties (Access = private)
        gidPath
        tclFileName
        filesFolder
        tclFileWriter
        tclTemplateName
        fieldName
        outputName
    end
    
    properties (Access = private)
        outPutNames
        fileName
        designVariableName     
        tclTemplateNames
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
            obj.tclTemplateNames   = cParams.tclTemplateNames;
            obj.outPutNames        = cParams.outPutNames;
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
            %obj.gidPath = '/home/alex/GiDx64/14.1.9d';
            obj.gidPath = '/home/alex/GiDx64/gid-15.0.3';
            obj.filesFolder = fullfile(pwd,'Output',obj.fileName);
        end
        
        function makeDesignVariableVideo(obj)
            for ivideo = 1:numel(obj.tclTemplateNames)
                obj.tclTemplateName = obj.tclTemplateNames{ivideo};
                obj.fieldName = obj.designVariableName;
                obj.outputName = [obj.fileName,obj.outPutNames{ivideo}];
                obj.makeFieldVideo();
            end
        end
        
        function makeRegularizedDesignVariableVideo(obj)
            obj.fieldName = 'DensityGauss';
            obj.tclTemplateName = 'Make_Video_density';
            obj.outputName = [obj.fileName,obj.fieldName];
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
            s.type            = obj.fieldName;
            s.tclFileName     = obj.tclFileName;
            s.filesFolder     = obj.filesFolder;
            s.outputName      = obj.outputName;
            s.iterations      = obj.iterations;
            s.fileName        = obj.fileName;
            s.tclTemplateName = obj.tclTemplateName;
            obj.tclFileWriter = TclFileWriter(s);
        end
        
        function writeTclFile(obj)
            obj.tclFileWriter.write();
        end
        
        function executeTclFile(obj)
            tFile = replace(obj.tclFileName,'\','\\');
            gFile = fullfile(obj.gidPath,'gid_offscreen');
            executingLine = ['"',gFile,'"',' -offscreen', ' -t ' ,'"source ',tFile,'"'];
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
