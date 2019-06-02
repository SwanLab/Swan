classdef VideoMaker_TopOpt < VideoMaker_Physical_Problem
    
    properties (Access = protected)
        fieldName  
    end
    
    properties (GetAccess = protected, SetAccess = private)
  
    end    
    
    methods (Access = public)
        
        function makeDesignVariableVideo(obj)
            obj.createVideoFileName();
            obj.createFinalPhotoName();
            obj.createFileList();          
            obj.createTclFileName();            
            obj.writeTclFileName();
            obj.executeTclFiles();    
            obj.deleteTclFile();
        end
        
        function makeRegDesignVariableVideo(obj)
            post = Postprocess_TopOpt_density();
            field2print = post.density_name_reg;
            componentfield = post.density_name_component_reg;
            obj.Make_video_density_field(field2print,componentfield,output_video_name)
        end          
        
    end       
    
    methods (Access = protected)
        
   
                
 
        
    end
    
    methods (Access = private)
        

    end

    
    methods (Access = protected, Abstract)       
      writeTclFileName(obj)        
    end    
    
end