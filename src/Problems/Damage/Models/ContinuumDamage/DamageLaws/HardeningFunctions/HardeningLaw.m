classdef HardeningLaw < handle
       
    properties (Access = protected)
        r0
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
            obj.r0 = cParams.r0;
        end

    end

end