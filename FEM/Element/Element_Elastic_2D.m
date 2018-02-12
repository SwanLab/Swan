classdef Element_Elastic_2D<Element_Elastic
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function variables = computeVars(obj,uL)
            variables = obj.computeDispStressStrain(uL);
            variables.strain = obj.computeEz(variables.strain,obj.nstre,obj.nelem,obj.material);
            variables = obj.permuteStressStrain(variables);
        end
    end
    
end

