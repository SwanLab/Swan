classdef QuadrilateralNodesCalculator < NodesCalculator

    properties (Access = public)
        nvert
        div
        boundNodes
        totalNodes
    end
    
    methods (Access = public)
        
        function obj = QuadrilateralNodesCalculator(cParams)
            obj.init(cParams);
            obj.computeBoundaryNodes();
            obj.computeTotalNodes();
        end
        
    end
    
    methods (Access = private)
        
        function computeTotalNodes(obj)
            divA = obj.div(1);
            divB = obj.div(2);
            obj.totalNodes = obj.boundNodes+(divA-1)*(divB-1);
        end
        
    end

end