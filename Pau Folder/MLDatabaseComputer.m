classdef MLDatabaseComputer < handle
    
    properties (Access = public)
        Chomog
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MLDatabaseComputer()
            radius = linspace(0,0.5,10);
            Chomog = zeros(3,3,length(radius));
            for i=1:length(radius)
                fem = Tutorial02p2FEMElasticityMicro(radius(i));
                figure()
                fem.mesh.plot
                obj.Chomog(:,:,i) = fem.stateProblem.Chomog;
            end

        end
        
    end
    
end