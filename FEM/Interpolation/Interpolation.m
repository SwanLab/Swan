classdef Interpolation < handle
    
    properties (GetAccess = public, SetAccess = protected)
        order
        
        ndime
        nnode
        
        pos_nodes
        % isoDv
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
        computeShapeFunctions(obj)
        computeShapeDerivatives(obj)
    end
    
end
