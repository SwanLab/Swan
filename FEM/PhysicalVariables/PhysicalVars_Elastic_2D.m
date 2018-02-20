classdef PhysicalVars_Elastic_2D < PhysicalVars_Elastic
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties
    end
    
    methods (Access = {?Physical_Problem, ?PhysicalVars_Elastic_2D, ?PhysicalVars_Elastic_2D_Micro})
        function obj = computeVars(obj,uL)
            obj.d_u(obj.dof.vL) = uL;
            obj.d_u(obj.dof.vR) = obj.bc.fixnodes;
            strain = obj.computeStrain(obj.d_u,obj.dim,obj.nnode,obj.nelem,obj.geometry.ngaus,obj.bc.idx);
            stress = obj.computeStress(strain,obj.material.C,obj.geometry.ngaus,obj.nstre);
            strain = obj.computeEz(strain,obj.nstre,obj.nelem,obj.material);
            obj.strain = strain;
            obj.stress = stress;            

            obj.strain = permute(strain, [3 1 2]);
            obj.stress = permute(stress, [3 1 2]);           
        end
    end
    
    methods (Access = protected, Static)
        % Compute strains

        function strain = computeEz(strain,dim,nelem,material)
            mu = full(material.mu);
            kappa = full(material.kappa);
            epoiss = (kappa(1,1) - mu(1,1))./(kappa(1,1) + mu(1,1));
            epoiss = ones(1,nelem)*epoiss;
            strain(nstre+1,:,:) = (-epoiss./(1-epoiss)).*(strain(1,:,:)+strain(2,:,:));
        end
    end
    
end

