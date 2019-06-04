classdef TclFileWriter < handle
    
    properties (Access = protected)        
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
    
    methods (Access = public, Static)
   
        function obj = create(cParams)
            f = TclFileWriterFactory();
            obj = f.create(cParams);
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