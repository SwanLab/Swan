
classdef PhysicalVars_Elastic_2D_Micro < PhysicalVars_Elastic_2D
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties (Access = {?Physical_Problem_Micro, ?PhysicalVars_Elastic_2D})
        stress_homog
        stress_fluct
        strain_fluct
    end
    
    methods (Access = {?Physical_Problem, ?PhysicalVars_Elastic_2D, ?PhysicalVars_Elastic_2D_Micro})
        function obj = PhysicalVars_Elastic_2D_Micro(ndof)
            obj.d_u = zeros(ndof,1);
        end
        
        function obj = computeVars(obj,d_u,dim,G,nelem,idx,element,material,vstrain)
            obj = computeVars@PhysicalVars_Elastic_2D(obj,d_u,dim,G,nelem,idx,element,material);
            
            obj.stress_fluct = obj.stress;
            obj.strain_fluct = obj.strain;
            Cmat = material.C;
            
            obj.stress = zeros(G.ngaus,dim.nstre,nelem);
            obj.strain = zeros(G.ngaus,dim.nstre,nelem);           
            obj.stress_homog = zeros(dim.nstre,1);
            vol_dom = sum(sum(G.dvolu)); 
            
            for igaus=1:G.ngaus
                obj.strain(igaus,1:dim.nstre,:) = vstrain.*ones(1,dim.nstre,nelem) + obj.strain_fluct(igaus,1:dim.nstre,:);
                for istre=1:dim.nstre
                    for jstre=1:dim.nstre
                        obj.stress(igaus,istre,:) = squeeze(obj.stress(igaus,istre,:)) + 1/vol_dom*squeeze(squeeze(Cmat(istre,jstre,:))).* squeeze(obj.strain(igaus,jstre,:));
                    end
                end
                % contribucion a la C homogeneizada
                for istre=1:dim.nstre
                    obj.stress_homog(istre) = obj.stress_homog(istre) +  1/vol_dom *(squeeze(obj.stress(igaus,istre,:)))'*G.dvolu(:,igaus);
                end                
            end
        end
    end
end

