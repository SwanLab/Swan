classdef ComposedGradFunction < L2Function
      
    properties (Access = private)
        f1
        f2
        operation
        quad
    end
    
    properties (Access = public)
        ndimf
    end
    
    methods (Access = public)
        
        function obj = ComposedGradFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xG)
            f1F = obj.f1.computeGradient(obj.quad);
            f2F = obj.f2.computeGradient(obj.quad);
            f1V = f1F.fValues;
            f2V = f2F.fValues;
            fxV = zeros(size(f1V,2),size(f1V,3));
            for i = 1:size(f1V,1)
                f1i = squeeze(f1V(i,:,:,:));
                f2i = squeeze(f2V(i,:,:,:));
                fxV = fxV + obj.operation(f1i,f2i);
            end
         end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.f1        = cParams.f1;
            obj.f2        = cParams.f2;
            obj.operation = cParams.operation;
            obj.ndimf     = cParams.ndimf;
            obj.quad      = cParams.quad;
        end
        
    end
    
end