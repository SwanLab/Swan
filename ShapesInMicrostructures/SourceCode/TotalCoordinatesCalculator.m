classdef TotalCoordinatesCalculator < handle
    
    properties (Access = public, Abstract)
        c
        theta
        nodes
        vertCoord
        boundCoord
        div
        totalCoord
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            obj = TotalCoordinatesCalculatorFactory.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.c  = cParams.c;
            obj.theta = cParams.theta;
            obj.nodes = cParams.nodes;
            obj.vertCoord = cParams.vertCoord;
            obj.boundCoord = cParams.boundCoord;
            obj.div = cParams.div;
            obj.totalCoord = [];
        end
        
        function initBoundary(obj)
            boundNodes = obj.nodes.bound;
            obj.totalCoord(1:boundNodes,:) = obj.boundCoord;
        end
        
    end
    
end