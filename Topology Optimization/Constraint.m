classdef Constraint < CC
    
    properties (Access = public)
        lambda
    end
    
    methods (Access = public)
        
        function obj = Constraint(settings,designVar,homogVarComputer,targetParameters)
            obj.init(settings, settings.constraint,designVar,homogVarComputer,targetParameters);
        end
        
        function updateFields(obj,iSF)
            obj.value(iSF,1)    = obj.shapeFunctions{iSF}.value;
            obj.gradient(:,iSF) = obj.shapeFunctions{iSF}.gradient;
        end
        
    end
    
end