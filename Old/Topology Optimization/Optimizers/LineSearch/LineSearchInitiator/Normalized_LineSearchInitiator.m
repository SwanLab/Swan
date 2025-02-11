classdef Normalized_LineSearchInitiator < LineSearchInitiator
    
    properties (Access = private)
        scalarProduct
    end
    
    methods (Access = public)
        
        function obj = Normalized_LineSearchInitiator(cParams)
            obj.init(cParams);
            obj.scalarProduct = cParams.scalarProduct;
        end
        
        function initStep = compute(obj,value)
            x = obj.designVariable.value;
            g = obj.objectiveFunction.gradient;
            xNorm = obj.scalarProduct.computeSP(x,x);
            gNorm = obj.scalarProduct.computeSP(g,g);
            initStep = 0.01*sqrt(xNorm/gNorm);
        end
        
    end
    
end

