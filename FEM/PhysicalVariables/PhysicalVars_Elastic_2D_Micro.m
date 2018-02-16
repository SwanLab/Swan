
classdef PhysicalVars_Elastic_2D_Micro < PhysicalVars_Elastic_2D
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties (Access = {?Physical_Problem_Micro, ?PhysicalVars_Elastic_2D})
        stress_homog
        stress_fluct
        strain_fluct
    end
    
    properties (GetAccess = public, SetAccess = ?Physical_Problem_Micro)
        Chomog
        tstress
        tstrain
    end
    
    methods (Access = {?Physical_Problem, ?PhysicalVars_Elastic_2D, ?PhysicalVars_Elastic_2D_Micro})
        function obj = PhysicalVars_Elastic_2D_Micro(ndof)
            obj.d_u = zeros(ndof,1);
        end
        
        function variables = computeVars(variables,d_u,dim,G,nelem,idx,element,material,vstrain)
            variables = computeVars@PhysicalVars_Elastic_2D(variables,d_u,dim,G,nelem,idx,element,material);
            
            variables.stress_fluct = variables.stress;
            variables.strain_fluct = variables.strain;
            Cmat = material.C;
            
            variables.stress = zeros(G.ngaus,dim.nstre,nelem);
            variables.strain = zeros(G.ngaus,dim.nstre,nelem);           
            variables.stress_homog = zeros(dim.nstre,1);
            vol_dom = sum(sum(G.dvolu)); 
            
            for igaus=1:G.ngaus
                variables.strain(igaus,1:dim.nstre,:) = vstrain.*ones(1,dim.nstre,nelem) + variables.strain_fluct(igaus,1:dim.nstre,:);
                for istre=1:dim.nstre
                    for jstre=1:dim.nstre
                        variables.stress(igaus,istre,:) = squeeze(variables.stress(igaus,istre,:)) + 1/vol_dom*squeeze(squeeze(Cmat(istre,jstre,:))).* squeeze(variables.strain(igaus,jstre,:));
                    end
                end
                % contribucion a la C homogeneizada
                for istre=1:dim.nstre
                    variables.stress_homog(istre) = variables.stress_homog(istre) +  1/vol_dom *(squeeze(variables.stress(igaus,istre,:)))'*G.dvolu(:,igaus);
                end                
            end
        end
    end
end

