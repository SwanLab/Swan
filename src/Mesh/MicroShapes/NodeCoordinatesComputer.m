classdef NodeCoordinatesComputer < handle
    
    properties (Access = private)
        c
        theta
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
            obj.c = cParams.c;
            obj.theta = cParams.theta;
            obj.div = cParams.div;
            obj.nodes = cParams.nodes;
        end
        
        function computeVertCoord(obj)
            s.c = obj.c;
            s.theta = obj.theta;
            s.nodes = obj.nodes;
            a = VertexCoordinatesCalculator(s);
            obj.vertCoord = a.vertCoord;
        end
        
        function computeBoundCoord(obj)
            s.c = obj.c;
            s.theta = obj.theta;
            s.nodes = obj.nodes;
            s.div = obj.div;
            s.vertCoord = obj.vertCoord;
            a = BoundaryCoordinatesCalculator(s);
            obj.boundCoord = a.boundCoord;
        end
        
        function computeTotalCoord(obj)
            s.c = obj.c;
            s.theta = obj.theta;
            s.nodes = obj.nodes;
            s.div = obj.div;
            s.vertCoord = obj.vertCoord;
            s.boundCoord = obj.boundCoord;
            a = TotalCoordinatesCalculator.create(s);
            obj.totalCoord = a.totalCoord;
        end
        
    end
    
    methods (Static)
        
        function pos = computeThePosition(c0,c,theta)
            pos = c0+c.*[cosd(theta) sind(theta)];
        end
        
    end

end