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
        shapeFun, shapeDer, shapeSecondDer
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

        function secondDeriv = computeShapeSecondDerivatives(obj,xV)
            secondDeriv = obj.shapeSecondDer.evaluate(xV);
        end

    end

    methods (Access = public)
        
        function init(obj,cParams)
            obj.type  = cParams.type;
            obj.order = cParams.order;
            obj.computeParams();
            obj.createShapeFunctions();
            obj.createShapeDerivatives();
            obj.createShapeSecondDerivatives();
        end

        function createShapeFunctions(obj)
            s.operation = @(xV) obj.evaluateShapeFunctions(xV);
            s.ndimf = 1;
            shape = DomainFunction(s);
            obj.shapeFun = shape;
        end

        function createShapeDerivatives(obj)
            s.operation = @(xV) obj.evaluateShapeDerivatives(xV);
            s.ndimf = 1;
            deriv = DomainFunction(s);
            obj.shapeDer = deriv;
        end

        function createShapeSecondDerivatives(obj)
            s.operation = @(xV) obj.evaluateShapeSecondDerivatives(xV);
            s.ndimf = 1;
            secondDeriv = DomainFunction(s);
            obj.shapeSecondDer = secondDeriv;
        end
        
    end
    
    methods (Abstract, Access = protected)
        computeParams(obj)
        evaluateShapeFunctions(obj)
        evaluateShapeDerivatives(obj)
        %evaluateShapeSecondDerivatives(obj) %% HAS TO BE IMPLEMENTED IN
        %ALL CLASSES
    end
    
end
