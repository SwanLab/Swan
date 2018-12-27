classdef LevelSetPrinter < Printer
    
    properties
    end
    
    methods (Access = public)
        
       function obj = LevelSetPrinter(quad,mesh)
            obj.init(quad,mesh)            
        end
        
    end
    
    methods (Access = private)
        
        function createPostProcess(obj)
            SetOpt = 'SLERP';
            obj.PostProcess = Postprocess_TopOpt.Create(SetOpt);
        end

    end
    
    
end

