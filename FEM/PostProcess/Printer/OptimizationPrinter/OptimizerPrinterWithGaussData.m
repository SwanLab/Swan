classdef OptimizerPrinterWithGaussData < OptimizerPrinter
    
    properties (Access = private)
        quad
    end
    
    
    methods (Access = protected)
        
        function obtainHasGaussDataAndQuad(obj)
            [h,q] = obj.obtainQuadHasGaussData(obj.cost.ShapeFuncs{1});
            obj.hasGaussData = h;
            obj.quad = q;
        end
        
        function createDataInputForCreateDataBase(obj,mesh,fileName)
            obj.createDataInputForCreateDataBase@OptimizerPrinter(mesh,fileName);
            if obj.hasGaussData
               obj.dI.quad = obj.quad;
            end
        end
        
    end
    
    methods (Access = private)
        
        function [h,q] = obtainQuadHasGaussData(obj,sh)
            h = true;
            phyPr = sh.getPhysicalProblem();
            q = phyPr.element.quadrature;
        end
    end
    
    
end