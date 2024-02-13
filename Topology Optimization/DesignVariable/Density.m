classdef Density < DesignVariable
    
    methods (Access = public)
        
        function obj = Density(cParams)
            obj.nVariables = 1;
            obj.init(cParams);
        end

        function fun = obtainDomainFunction(obj)
            fun = obj.fun;
        end        

        function update(obj,value)
            if ~isempty(obj.isFixed)
                value(obj.isFixed.nodes) = obj.isFixed.values;
            end
            s.mesh    = obj.mesh;
            s.fValues = value;
            s.order   = 'P1';
            obj.fun   = LagrangianFunction(s);
        end        
    
    end
    
end

