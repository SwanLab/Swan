classdef PhysicalVars_Elastic_2D < PhysicalVars_Elastic
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function obj = computeVars(obj,d_u,dim,nnode,nelem,ngaus,idx,element,material)
            nstre = 3;
            obj.d_u = d_u;
            strain = obj.computeStrain(d_u,dim,nnode,nelem,ngaus,idx,element,material);
            obj.strain = strain;
            obj.stress = obj.computeStress(strain,material.C,ngaus,nstre);
        end        
    end
    
    methods (Access = protected, Static)
        % Compute strains
        function strain = computeStrain(d_u,dim,nnode,nelem,ngaus,idx,element,material)
            strain = computeStrain@PhysicalVars_Elastic(d_u,dim,nnode,nelem,ngaus,idx,element);
            
            % Compute eZ
            mu = material.mu;
            kappa = material.kappa;
            epoiss = (kappa - mu)./(kappa + mu);
            epoiss = ones(1,nelem)*epoiss;
            strain(dim.nstre+1,:,:) = (-epoiss./(1-epoiss)).*(strain(1,:,:)+strain(2,:,:));
        end
        
        % Compute stresses
        function stress = computeStress(strain,C,ngaus,nstre)
            stress = computeStress@PhysicalVars_Elastic(strain,C,ngaus,nstre);
        end
    end
    
end

