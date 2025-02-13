classdef TotalCoordinatesCalculatorFactory < handle
    methods (Access = public, Static) 
        function obj = create(cParams)
            switch cParams.nodes.vert
                case 4
                    obj = IntersectionCoordComputer(cParams);
                case 6
                    obj = DiagonalCoordComputer(cParams);
            end
        end
    end 
end