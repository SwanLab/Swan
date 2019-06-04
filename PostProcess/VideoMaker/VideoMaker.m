classdef VideoMaker < handle
    
    properties (GetAccess = private, SetAccess = public)
        iterations        
    end
    
    properties (GetAccess = protected, SetAccess = private)
        videoFileName
        fileList
        photoFileName         
    end
    
    properties (Access = protected)
        gidPath
        tclFileName
        filesFolder
        fileName
        tclFileWriter
        fieldName 
    end
    
    methods (Access = public)
        
        function makeVideo(obj,nIter)
            obj.iterations = 0:nIter;
            obj.makeDesignVariableVideo();
            obj.makeRegDesignVariableVideo(); 
        end
        
    end
    
    methods (Access = public, Static)
        
        function obj = VideoMaker(cParams)
            obj.init(cParams);
        end
        
        function obj = create(cParams)
            f = VideoMakerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.fileName    = cParams.caseFileName;
            obj.fieldName   = cParams.designVarType;
            obj.createPaths(cParams);
            obj.createFolder();              
        end
        
        function createVideoFileName(obj)
            iterStr = int2str(obj.iterations(end));
            fName = ['Video_',obj.fileName,'_',iterStr,'.gif'];
            fullName = fullfile(pwd,'Output',obj.fileName,fName);
            obj.videoFileName = SpecialCharacterReplacer.replace(fullName);
        end
        
        function createFileList(obj)
            iter2print = obj.iterations;
            fName      = obj.fileName;
            folderName = obj.filesFolder;
            list = [];
            for iter = 1:length(iter2print)
                iStr          = num2str(iter2print(iter));
                iFileName     = [fName,iStr,'.flavia.res'];
                iFullFileName = fullfile(folderName,iFileName);
                list = [list, ' ',iFullFileName];
            end
            obj.fileList = SpecialCharacterReplacer.replace(list);
        end
        
        function createFinalPhotoName(obj)
            fName = fullfile(obj.filesFolder,[obj.fileName,'.png']);
            obj.photoFileName = SpecialCharacterReplacer.replace(fName);
        end
        
        function createTclFileName(obj)
            fName = 'tcl_gid.tcl';
            obj.tclFileName = fullfile(obj.filesFolder,fName);
        end          
        
        function createTclFileWriter(obj,field)
            cParams.type          = field;
            cParams.tclFileName   = obj.tclFileName;            
            cParams.fileList      = obj.fileList;
            cParams.videoFileName = obj.videoFileName;
            cParams.photoFileName = obj.photoFileName;
            obj.tclFileWriter = TclFileWriter.create(cParams);
        end
        
        function writeTclFile(obj)
           obj.tclFileWriter.write(); 
        end
        
        function executeTclFiles(obj) 
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
    
    methods (Access = private)
        
       function createFolder(obj)
            if ~exist(obj.filesFolder,'dir')
                mkdir(obj.filesFolder);
            end            
       end     
        
        function createPaths(obj,cParams)
            fName = cParams.caseFileName;
            %obj.gidPath = 'C:\Program Files\GiD\GiD 13.0.4';% 
            obj.gidPath = '/opt/GiDx64/13.0.2/';
            obj.fileName = fName;
            obj.filesFolder = fullfile(pwd,'Output',fName);
        end       
        
        function makeDesignVariableVideo(obj)
            obj.createVideoFileName();
            obj.createFinalPhotoName();
            obj.createFileList();          
            obj.createTclFileName();  
            obj.createTclFileWriter(obj.fieldName);
            obj.writeTclFile();
            obj.executeTclFiles();    
            obj.deleteTclFile();
        end
        
        function makeRegDesignVariableVideo(obj)
            obj.createVideoFileName();
            obj.createFinalPhotoName();
            obj.createFileList();          
            obj.createTclFileName();  
            obj.createTclFileWriter('RegularizedDensity');
            obj.writeTclFile();
            obj.executeTclFiles();    
            obj.deleteTclFile();
        end               
        
        
        
    end


end
