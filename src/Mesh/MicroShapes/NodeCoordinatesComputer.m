classdef NodeCoordinatesComputer < handle
    
    properties (Access = private)
        latticeVectors
        div
        nodes
    end
    
    properties (Access = public)
        vertCoord
        boundCoord
        totalCoord
    end
    
    methods (Access = public)
        
        function obj = NodeCoordinatesComputer(cParams)
            obj.init(cParams); 
        end
        
        function computeCoordinates(obj)
            obj.computeVertCoord();
            obj.computeBoundCoord();
            obj.computeTotalCoord();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.latticeVectors = cParams.latticeVectors;
            obj.div = cParams.div;
            obj.nodes = cParams.nodes;
        end
        
        function computeVertCoord(obj)
            s.latticeVectors = obj.latticeVectors;
            s.nodes = obj.nodes;
            a = VertexCoordinatesCalculator(s);
            obj.vertCoord = a.vertCoord;
        end
        
        function computeBoundCoord(obj)
            s.latticeVectors = obj.latticeVectors;
            s.nodes = obj.nodes;
            s.div = obj.div;
            s.vertCoord = obj.vertCoord;
            a = BoundaryCoordinatesCalculator(s);
            obj.boundCoord = a.boundCoord;
        end
        
        function computeTotalCoord(obj)
            s.latticeVectors = obj.latticeVectors;
            s.nodes = obj.nodes;
            s.div = obj.div;
            s.vertCoord = obj.vertCoord;
            s.boundCoord = obj.boundCoord;
            a = TotalCoordinatesCalculator.create(s);
            obj.totalCoord = a.totalCoord;
        end
        
    end
    
end