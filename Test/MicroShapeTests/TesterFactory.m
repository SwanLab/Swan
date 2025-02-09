classdef TesterFactory < handle
    methods (Access = public, Static) 
        function obj = create(type,initialData)
            switch type
                case 'NodeCoordinatesComputerTester'
                    obj = NodeCoordinatesComputerTester(initialData);
                case 'MeshCreatorTester'
                    obj = MeshCreatorTester(initialData);
            end
        end
    end
end