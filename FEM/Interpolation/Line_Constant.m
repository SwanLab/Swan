classdef Line_Constant < Interpolation

    properties (Access = private)
        ngaus
    end

    methods (Access = public)

        function obj = Line_Constant(cParams)
            obj.init(cParams);
            obj.computeParams();
        end

        function computeShapeDeriv(obj,posgp)
            obj.ngaus = length(posgp);
            obj.computeShapes()
            obj.computeShapeDerivatives();
        end

    end

    methods (Access = private)

        function computeParams(obj)
            obj.type = 'LINE';
            obj.ndime = 2;
            obj.nnode = 1;
            obj.pos_nodes = [0 0];
        end

        function computeShapes(obj)
            N = ones(1,obj.ngaus);
            obj.shape = N;
        end

        function computeShapeDerivatives(obj)
            dN = zeros(2,obj.ngaus);
            obj.deriv = dN;
        end
        
    end
    
end