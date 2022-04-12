classdef TesterFactory < handle
    methods (Access = public, Static) 
        function obj = create(type,initialData)
            switch type
                case 'NodesCalculatorTester'
                    obj = NodesCalculatorTester(initialData);
                case 'QuadrilateralNodesCalculatorTester'
                    obj = QuadrilateralNodesCalculatorTester(initialData);
                case 'HexagonalNodesCalculatorTester'
                     obj = HexagonalNodesCalculatorTester(initialData);
                case 'VertexCoordinatesCalculatorTester'
                    obj = VertexCoordinatesCalculatorTester(initialData);
                case 'BoundaryCoordinatesCalculatorTester'
                    obj = BoundaryCoordinatesCalculatorTester(initialData);
                case 'NodeCoordinatesComputerTester'
                    obj = NodeCoordinatesComputerTester(initialData);
                case 'IntersectionCoordComputerTester'
                    obj = IntersectionCoordComputerTester(initialData);
                case 'DiagonalCoordComputerTester'
                     obj = DiagonalCoordComputerTester(initialData);
                case 'MasterSlaveComputerTester'
                    obj = MasterSlaveComputerTester(initialData);
                case 'MeshCreatorTester'
                    obj = MeshCreatorTester(initialData);
            end
        end
    end
end