classdef Density < DesignVariable
    
    methods (Access = public)
        
        function obj = Density(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.fun.fValues;
        end
                       
    end
    
end

