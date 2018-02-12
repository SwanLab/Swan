classdef Element_Elastic_3D<Element_Elastic
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function variables = computeVars(obj,uL)
            variables = obj.computeDispStressStrain(uL);
            variables = obj.permuteStressStrain(variables);
        end
    end
    
end



