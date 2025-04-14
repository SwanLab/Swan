classdef TrainedRVE < handle
    
    properties (Access = public)
        ndimf
        Kcoarse
        Udef
        Urb
        PhiDef
        PhiRb
        Grb
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
            obj.PhiDef  = EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS;
            obj.PhiRb   = EIFEoper.RECONSTRUCTION.RB_DISP.BASIS;
%             obj.Grb     = EIFEoper.Grb;
        end
        
    end
    
end