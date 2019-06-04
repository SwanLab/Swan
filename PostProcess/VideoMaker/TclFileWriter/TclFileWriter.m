classdef TclFileWriter < handle
    
    properties (Access = protected)        
      fieldName
      tclFileName      
      tclTemplateName      
      fileList
      videoFileName
      photoFileName

      fullTclTemplateName     
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
           obj.fileList        = cParams.fileList;
           obj.videoFileName   = cParams.videoFileName;
           obj.photoFileName   = cParams.photoFileName;
           obj.tclFileName     = cParams.tclFileName;
           obj.createFullTclTemplateName();
        end
        
    end
        
    
    methods (Access = private)
        
       function createFullTclTemplateName(obj)
            fName = fullfile(pwd,'PostProcess','VideoMaker',[obj.tclTemplateName,'.tcl']);
            fName = SpecialCharacterReplacer.replace(fName);
            obj.fullTclTemplateName = fName;
       end        
        
       
    
    end
    
    
end