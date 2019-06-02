classdef VideoMaker < handle
    
    properties (GetAccess = protected, SetAccess = private)
        videoFileName
        fullTclTemplateName     
        fileList
        photoFileName         
    end
    
    properties (Access = protected)
        gidPath
        tclFileName
        tclTemplateName                  
        filesFolder
        iterations
        fileName
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = VideoMakerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function Set_up_make_video(obj,iter)
            obj.iterations = iter;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.gidPath     = cParams.gidPath;
            obj.fileName    = cParams.fileName;
            obj.filesFolder = cParams.filesFolder;
            obj.createTclTemplateName();
        end
        
        function createVideoFileName(obj)
            iterStr = int2str(obj.iterations(end));
            fName = ['Video_',obj.fileName,'_',iterStr,'.gif'];
            fullName = fullfile(pwd,'Output',obj.fileName,fName);
            obj.videoFileName = obj.replace_special_character(fullName);
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
            obj.fileList = obj.replace_special_character(list);
        end
        
        function createFinalPhotoName(obj)
            fName = fullfile(obj.filesFolder,[obj.fileName,'.png']);
            obj.photoFileName = obj.replace_special_character(fName);
        end
        
        function createTclFileName(obj)
            fName = 'tcl_gid.tcl';
            obj.tclFileName = fullfile(obj.filesFolder,fName);
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
    
    methods (Access = protected, Static)
        
        function [output_string] = replace_special_character(input_string)
            if ispc
                output_string = replace(input_string,'\','\\\\');
            elseif isunix
                output_string = input_string;
            elseif ismac
            end
        end
        
    end
    
    methods (Access = private)
        
       function createTclTemplateName(obj)
            fName = fullfile(pwd,'PostProcess','VideoMaker',[obj.tclTemplateName,'.tcl']);
            obj.fullTclTemplateName = obj.replace_special_character(fName);
        end           
        
      
        
    end
end
