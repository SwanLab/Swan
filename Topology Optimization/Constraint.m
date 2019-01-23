classdef Constraint < CC
    properties (Access = public)
        lambda
    end
    
    methods (Access = public)
        function obj=Constraint(settings)
            obj@CC(settings,settings.constraint);
        end
        
        function updateFields(obj,iSF)
            obj.value(iSF,1) = obj.ShapeFuncs{iSF}.value;
            obj.gradient(:,iSF) = obj.ShapeFuncs{iSF}.gradient;
        end
    end
end