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
    
    methods (Access = protected, Static)
        
        function results = createResultsInputStructure(ls,outname)
            results.iter = 0;
            results.case_file = outname;
            results.design_variable = ls;
        end
        
    end
    
end

