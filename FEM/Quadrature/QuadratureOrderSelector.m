classdef QuadratureOrderSelector < handle

    properties (Access = public)
        
    end

    properties (Access = private)
        func
        ndim
    end

    methods (Access = public)

        function obj = QuadratureOrderSelector(cParams)
            obj.func = cParams.func;
            obj.ndim = cParams.ndim;
        end
        
    end
    
    methods (Access = private)

        function quadratureOrder = computeOrderLagrangian(obj,func)
            if func.order(1) == 'P'
                funcOrder = str2double(func.order(2));
            else
                funcOrder = str2double(func.order(6));
            end

            quadratureOrder = ceil(funcOrder/obj.ndim);
        end

        function quadratureOrder = orderLiteral(~,ord)
            switch ord
                case 0
                    quadratureOrder = 'CONSTANT';
                case 1
                    quadratureOrder = 'LINEAR';
                case 2
                    quadratureOrder = 'QUADRATIC';
                case 3
                    quadratureOrder = 'CUBIC';
                otherwise
                    quadratureOrder = 'ORDER';
                    quadratureOrder(end+1) = num2str(ord);
            end
        end

    end

    methods (Access = public, Static)

        function quadratureOrder = compute(cParams)

            obj = QuadratureOrderSelector(cParams);

            switch class(obj.func)
                case 'LagrangianFunction'
                    quadratureOrder = obj.computeOrderLagrangian(obj.func);
            end

            quadratureOrder = obj.orderLiteral(quadratureOrder);

        end

    end

end