classdef Shape_Functional < handle
    properties
        value
        gradient
        target_parameters=struct;
    end    
    methods
        function obj=Shape_Functional(settings)
            obj.target_parameters=settings.target_parameters;            
        end
        computef(obj)
    end
end
