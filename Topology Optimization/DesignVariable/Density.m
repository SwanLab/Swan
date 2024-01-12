classdef Density < DesignVariable
    
    methods (Access = public)
        
        function obj = Density(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
        end

        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            s.mesh    = obj.mesh;
            s.fValues = value;
            obj.fun   = P1Function(s);
        end        
    
    end
    
end

