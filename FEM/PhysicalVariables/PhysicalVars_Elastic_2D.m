classdef PhysicalVars_Elastic_2D < PhysicalVars_Elastic
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties
    end
    
    methods (Access = {?Physical_Problem, ?PhysicalVars_Elastic_2D, ?PhysicalVars_Elastic_2D_Micro})
        function obj = computeVars(obj,d_u,dim,G,nelem,idx,element,material)
            obj.d_u = d_u;
            strain = obj.computeStrain(d_u,dim,G.nnode,nelem,G.ngaus,idx,element);
            stress = obj.computeStress(strain,material.C,G.ngaus,dim.nstre);
            strain = obj.computeEz(strain,dim,nelem,material);
           
            obj.strain = permute(strain, [3 1 2]);
            obj.stress = permute(stress, [3 1 2]);           
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

