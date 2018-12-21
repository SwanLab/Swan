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
    
    methods (Access = protected, Static)
        
        function results = createResultsInputStructure(dens,outname,iter)
            results.iter = iter;
            results.case_file = strcat(outname,'Density');
            results.density = dens;
        end
    end
    
end

