classdef TotalCoordinatesCalculator < handle
    
    properties (Access = public, Abstract)
        latticeVectors
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
            
            if isfield(cParams, 'latticeVectors')
                obj.latticeVectors = cParams.latticeVectors;
            end
            if isfield(cParams, 'c')
                obj.c = cParams.c;
            end
            if isfield(cParams, 'theta')
                obj.theta = cParams.theta;
            end
            obj.nodes = cParams.nodes;
            obj.vertCoord = cParams.vertCoord;
            obj.boundCoord = cParams.boundCoord;
            obj.div = cParams.div;
            obj.totalCoord = zeros(cParams.nodes.total, 2);
        end
        
        function initBoundary(obj)
            boundNodes = obj.nodes.bound;
            obj.totalCoord(1:boundNodes,:) = obj.boundCoord;
        end
        
    end
    
end