classdef Triangle_Linear < Interpolation

    properties (Access = private)
        shapeFun
        shapeDer
    end

    methods (Access = public)

        function obj = Triangle_Linear(cParams)
            obj.createParams();
            obj.createShapeFunctions();
            obj.createShapeDerivatives();
        end

        function shape = computeShapeFunctions(obj,xV)
            shape = obj.shapeFun.evaluate(xV);
        end

        function deriv = computeShapeDerivatives(obj,xV)
            deriv = obj.shapeDer.evaluate(xV);
        end

    end

    methods (Access = private)

        function createParams(obj)
            obj.ndime     = 2;
            obj.nnode     = 3;
            obj.pos_nodes = [0 0; 1 0; 0 1];
        end

        function createShapeFunctions(obj)
            s.operation = @(xV) obj.evaluateShapeFun(xV);
            shape = DomainFunction(s);
            obj.shapeFun = shape;
        end

        function createShapeDerivatives(obj)
            s.operation = @(xV) obj.evaluateShapeDer(xV);
            deriv = DomainFunction(s);
            obj.shapeDer = deriv;

        end

        function shape = evaluateShapeFun(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            s = xV(1,:,:);
            t = xV(2,:,:);
            I = ones(size(t));
            shape = zeros(obj.nnode,ngaus,nelem);
            shape(1,:,:) = I-s-t;
            shape(2,:,:) = s;
            shape(3,:,:) = t;
        end

        function deriv = evaluateShapeDer(obj, xV)
            ngaus = size(xV,2);
            nelem = size(xV,3);
            deriv = zeros(obj.ndime,obj.nnode,ngaus,nelem);
            deriv(1,1,:,:) = -1;
            deriv(1,2,:,:) = 1;
            deriv(1,3,:,:) = 0;
            deriv(2,1,:,:) = -1;
            deriv(2,2,:,:) = 0;
            deriv(2,3,:,:) = 1;
        end

    end

end
