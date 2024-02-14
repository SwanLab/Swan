classdef QuadratureOrderSelector < handle

    properties (Access = public)
        
    end

    properties (Access = private)
        funcs
        ndim
        mType
    end

    methods (Access = public)

        function obj = QuadratureOrderSelector(cParams)
            obj.funcs = cParams.funcs;
            obj.ndim = cParams.ndim;
            obj.mType = cParams.mType;
        end
        
    end
    
    methods (Access = private)

        function quadratureOrder = computeOrderLagrangian(obj,func)

                if func.order(1) == 'P'
                    funcOrder = str2double(func.order(2));
                else
                    funcOrder = str2double(func.order(6));
                end

                switch obj.mType
                    case {'TETRAHEDRA','TRIANGLE'}
                        quadratureOrder = funcOrder;
                    case {'QUAD','LINE','HEXAHEDRA'}
                        quadratureOrder = funcOrder*obj.ndim;
                end
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
                    quadratureOrder(end+1:end+length(num2str(ord))) = num2str(ord);
            end
        end

    end

    methods (Access = public, Static)

        function quadratureOrder = compute(cParams)

            obj = QuadratureOrderSelector(cParams);
            quadratureOrder = zeros(size(cParams.funcs));

            for iFunc = 1:length(obj.funcs)
                switch class(obj.funcs{iFunc})
                    case 'LagrangianFunction'
                        quadratureOrder(iFunc) = obj.computeOrderLagrangian(obj.funcs{iFunc});

                    case 'AnalyticalFunction'
                        quadratureOrder(iFunc) = 2;
    
                end
            end

            quadratureOrder = obj.orderLiteral(sum(quadratureOrder));

        end

    end

end