classdef Material_Elastic<Material
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = Material_Elastic()
            obj.matProp.kappa = .9107;
            obj.matProp.mu = .3446;
        end
    end
    
end

