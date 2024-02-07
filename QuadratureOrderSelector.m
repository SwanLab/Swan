classdef QuadratureOrderSelector < handle
    
    properties (Access = public)
        quadratureOrder
        ndim
        type
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = QuadratureOrderSelector(cParams)
            obj.init(cParams);
        end

        function quadratureOrder = computeOrder(obj,func)
            switch func.fType
                case 'FE'
                    quadratureOrder = obj.computeOrderLagrangian(func);
            end

            quadratureOrder = obj.orderLiteral(quadratureOrder);

        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.type = cParams.type;
            obj.ndim = cParams.ndim;
        end

        function quadratureOrder = computeOrderLagrangian(obj,func)
            funcOrder = str2double(func.order(2));
            quadratureOrder = funcOrder*obj.ndim;
        end

        function quadratureOrder = orderLiteral(~,ord)
            switch ord
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

end