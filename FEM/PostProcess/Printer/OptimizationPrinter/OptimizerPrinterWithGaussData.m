classdef OptimizerPrinterWithGaussData < OptimizerPrinter
    
    properties (Access = private)
        quad
    end
    
    
    methods (Access = protected)
        
        function obtainHasGaussDataAndQuad(obj)
            [h,q] = obj.obtainQuadHasGaussData(obj.cost.ShapeFuncs{1});
            obj.hasGaussData = h;
            obj.quad = q;
            if obj.hasGaussData
                obj.dI.quad = obj.quad;
            end
        end
        
    end
    
    
end