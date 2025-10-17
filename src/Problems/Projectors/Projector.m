classdef Projector < handle

  

    methods (Static, Access = public)

        function obj = create(cParams)
            p = ProjectorFactory();
            obj = p.create(cParams);
        end

    end

    methods (Static, Access = protected)

        function ord = determineQuadratureOrder(fun)
            switch class(fun)
                case 'ConstantFunction'
                    ord = 0;
                case 'AnalyticalFunction'
                    ord = 2;
                case 'FEFunction'
                    ord = 2; % needed to project P1 to P1D             
                otherwise
                    ord = 3;
            end
        end

    end

end