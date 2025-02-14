classdef TclFileWriter < handle
    
    properties (Access = private)        
      fieldName
      tclFileName      
      tclTemplateName      
      fileList
      filesFolder
      videoFileName
      photoFileName
      iterations
      fileName
      fullTclTemplateName     
      outputName      
    end
        
    methods (Access = public) 

        function obj = TclFileWriter(cParams)
           obj.init(cParams);           
        end        
        
        function write(obj)
            fid = fopen(obj.tclFileName,'w+');            
            fprintf(fid,'GiD_Process PostProcess \n');
            fprintf(fid,'%s\n',['set arg1 "',obj.fileList,'"']);
            fprintf(fid,'%s\n',['set arg2 "',obj.videoFileName,'"']);
            fprintf(fid,'%s\n',['set arg3 "',obj.fieldName,'"']);
            fprintf(fid,'%s\n',['set arg4 "',obj.fieldName,'"']);
            fprintf(fid,'%s\n',['set arg5 "',obj.photoFileName,'"']);            
            fprintf(fid,'%s\n',['source "',obj.fullTclTemplateName,'"']);
            fprintf(fid,'%s\n',[obj.tclTemplateName,' $arg1 $arg2 $arg3 $arg4 $arg5']);
            fprintf(fid,'%s\n','GiD_Process Mescape Quit');
            fclose(fid);         
        end         
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
           obj.fieldName       = cParams.type;
           obj.filesFolder     = cParams.filesFolder;
           obj.iterations      = cParams.iterations;
           obj.fileName        = cParams.fileName;
           obj.tclFileName     = cParams.tclFileName;
           obj.outputName      = cParams.outputName;
           obj.tclTemplateName = cParams.tclTemplateName;                                           
           obj.createVideoFileName();
           obj.createFileList();           
           obj.createFinalPhotoName();
           obj.createFullTclTemplateName();
        end
        
    end        
    
    methods (Access = private)
        
       function createFullTclTemplateName(obj)
            fName = fullfile(pwd,'PostProcess','VideoMaker',[obj.tclTemplateName,'.tcl']);
            fName = SpecialCharacterReplacer.replace(fName);
            obj.fullTclTemplateName = fName;
       end       
       
        function createVideoFileName(obj)
            iterStr = int2str(obj.iterations(end));
            fName = ['Video_',obj.outputName,'_',iterStr,'.gif'];
            fullName = fullfile(obj.filesFolder,fName);
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
            fName = fullfile(obj.filesFolder,[obj.outputName,'.png']);
            obj.photoFileName = SpecialCharacterReplacer.replace(fName);
        end       
        
    end
    
    
end