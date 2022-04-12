classdef NodesCalculatorFactory < handle
    methods (Access = public, Static) 
        function obj = create(cParams)
            switch cParams.nvert
                case 4
                    obj = QuadrilateralNodesCalculator(cParams);
                case 6
                    obj = HexagonalNodesCalculator(cParams);
            end
        end
    end
end