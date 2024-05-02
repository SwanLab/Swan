classdef TrainedRVE < handle
    
    properties (Access = public)
        ndimf
        Kcoarse
        Udef
        Urb
    end
    
    properties (Access = private)
        
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = TrainedRVE(filename)
            obj.init(filename)
            
        end
                
    end
    
    methods (Access = private)
        
        function init(obj,filename)
            load(filename);
            obj.Kcoarse = EIFEoper.Kcoarse;
            obj.Udef    = EIFEoper.Udef;
            obj.Urb     = EIFEoper.Urb;
            obj.ndimf   = 2;
        end
        
    end
    
end