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
            strain = obj.computeStrain(d_u,dim,nnode,nelem,ngaus,idx,element);
            strain = obj.computeEz(strain,dim,nelem,material);
            obj.strain = permute(strain, [3 1 2]);
            obj.stress = obj.computeStress(strain,material.C,ngaus,nstre);
            obj.stress = permute(obj.stress, [3 1 2]);
        end
    end
    
    methods (Access = protected, Static)
        % Compute strains
        function strain = computeEz(strain,dim,nelem,material)
            mu = material.mu;
            kappa = material.kappa;
            epoiss = (kappa(1,1) - mu(1,1))./(kappa(1,1) + mu(1,1));
            epoiss = ones(1,nelem)*epoiss;
            strain(dim.nstre+1,:,:) = (-epoiss./(1-epoiss)).*(strain(1,:,:)+strain(2,:,:));
        end
    end
    
end

