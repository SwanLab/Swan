classdef CoarseTrainedRVE < handle
    
    properties (Access = public)
        ndimf
        Kcoarse
        U
    end
    
    properties (Access = private)
        
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CoarseTrainedRVE(filename)
            obj.init(filename)          
        end
                
    end
    
    methods (Access = private)
        
        function init(obj,filename)
            data = load(filename);
            obj.Kcoarse = data.L;
            obj.U       = data.U;
            obj.ndimf   = 2;
        end
        
    end
    
end