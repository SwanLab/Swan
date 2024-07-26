classdef Display_SemilogY < Display_Abstract

    methods (Access = public)

        function obj = Display_SemilogY(cParams)
            obj@Display_Abstract(cParams.title)
        end

    end

    
    methods (Access = protected)
        
        function setChartType(obj)
            obj.handle = semilogy(0,0);
        end
        
    end
    
end