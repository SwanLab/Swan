
classdef PhysicalVars_Elastic_2D_Micro < PhysicalVars_Elastic_2D
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties
        stress_homog
    end
    
    methods (Access = ?Physical_Problem)
        function obj = computeVars(obj,d_u,dim,G,nelem,idx,element,material)
            computeVars@PhysicalVars_Elastic_2D(d_u,dim,G,nelem,idx,element,material);
            
            strain = obj.computeStrain(d_u,dim,G.nnode,nelem,G.ngaus,idx,element);
            strain = obj.computeEz(strain,dim,nelem,material);
            obj.strain = permute(strain, [3 1 2]);
            stress = obj.computeStress(strain,material.C,G.ngaus,dim.nstre);
            stress = permute(stress, [3 1 2]);
            
            for i=1:size(stress,2)
                newStress(i) = sum(stress(:,i,:))*G.dvolu;
            end
            obj.stress_homog = newStress./sum(G.dvolu);
            
        end
        
    end
end
