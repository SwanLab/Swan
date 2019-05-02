classdef Constraint < CC
    
    properties (Access = public)
        lambda
    end
    
    properties (Access = private)
       dualVariable 
    end
    
    methods (Access = public)
        
        function obj = Constraint(cParams)
            cParams.shapeFuncList = cParams.settings.constraint;
            obj.dualVariable = cParams.dualVariable;
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function updateFields(obj,iSF)
            obj.value(iSF,1)    = obj.shapeFunctions{iSF}.value;
            obj.gradient(:,iSF) = obj.shapeFunctions{iSF}.gradient;
        end
        
    end
    
end