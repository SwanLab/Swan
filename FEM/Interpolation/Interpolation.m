classdef Interpolation < handle
    
    properties (GetAccess = public, SetAccess = protected)
      
        order
        
        ndime
        nnode
        
        pos_nodes
        shape
        deriv
        isoDv
        
        iteration
        cases
        selectcases
        main_loop
        extra_cases
    end
    
    properties (Access = protected)
        type
    end
    
    methods (Static, Access = public)
        
        function obj = create(type,order)
            cParams.type = type;
            cParams.order = order;
            f = InterpolationFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type  = cParams.type;
            obj.order = cParams.order;
        end
        
    end
    
    methods (Abstract)
        computeShapeDeriv(obj)
    end
end
