classdef HardeningLaw < handle
       
    properties (Access = private)
        r
        r0
        isDamaging
        isOverR1
    end
    
    methods (Access = public, Static)
        function obj = create(s)
            f = HardeningLawFactory();
            obj = f.create(s); 
        end
    end
    
    methods (Access = public)
        
        function obj = HardeningLaw(cParams)
            obj.init(cParams);
        end
    end
    
    methods (Access = protected)
        function init(obj,cParams)
            % obj.r = cParams.r;
            % obj.r0 = cParams.r0;
            % obj.isDamaged = cParams.isDamaged;
            % obj.isOverR1 = cParams.isOverR1;

            obj.r = @(xV) cParams.r.evaluate(xV);
            obj.r0 = @(xV) cParams.r0.evaluate(xV);
            obj.isDamaging = @(xV) cParams.isDamaging.evaluate(xV);
            obj.isOverR1 = @(xV) cParams.isOverR1.evaluate(xV);
        end 
    end
end