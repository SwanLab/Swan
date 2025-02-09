classdef NodesCalculator < handle
    
    properties (Access = public, Abstract)
        nvert
        div
        boundNodes
        totalNodes
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            obj = NodesCalculatorFactory.create(cParams);
        end

    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.nvert = cParams.nvert;
            obj.div = cParams.div;
        end
        
        function computeBoundaryNodes(obj)
            halfIntNodes = 0;
            for iNode = 1:length(obj.div)
                halfIntNodes = halfIntNodes+obj.div(iNode)-1;
            end
            obj.boundNodes = round(obj.nvert+2*halfIntNodes);
        end
        
    end

end