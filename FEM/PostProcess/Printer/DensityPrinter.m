classdef DensityPrinter < Printer
    

    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = DensityPrinter(quad,mesh)
            obj.init(quad,mesh)            
        end

    end
    
    methods (Access = protected)
        
        function createPostProcess(obj)
            quad = obj.quadrature;
            obj.PostProcess = PostprocessDensityInGaussPoints(quad);
        end
       
    end
    
   
end

