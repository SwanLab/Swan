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
        shapeFun, shapeDer
    end
    
    methods (Static, Access = public)
        
        function obj = create(type,order)
            cParams.type = type;
            cParams.order = order;
            f = InterpolationFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)

        function shape = computeShapeFunctions(obj,xV)
            shape = obj.shapeFun.evaluate(xV);
        end

        function deriv = computeShapeDerivatives(obj,xV)
            deriv = obj.shapeDer.evaluate(xV);
        end

    end

    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type  = cParams.type;
            obj.order = cParams.order;
            obj.computeParams();
            obj.createShapeFunctions();
            obj.createShapeDerivatives();
        end

        function createShapeFunctions(obj)
            s.operation = @(xV) obj.evaluateShapeFunctions(xV);
            shape = DomainFunction(s);
            obj.shapeFun = shape;
        end

        function createShapeDerivatives(obj)
            s.operation = @(xV) obj.evaluateShapeDerivatives(xV);
            deriv = DomainFunction(s);
            obj.shapeDer = deriv;
        end
        
    end
    
    methods (Abstract, Access = protected)
        computeParams(obj)
        evaluateShapeFunctions(obj)
        evaluateShapeDerivatives(obj)
    end
    
end
