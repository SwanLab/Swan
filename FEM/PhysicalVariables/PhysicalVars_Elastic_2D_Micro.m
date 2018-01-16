
classdef PhysicalVars_Elastic_2D_Micro < PhysicalVars_Elastic_2D
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties (Access = {?Physical_Problem_Micro, ?PhysicalVars_Elastic_2D})
        stress_homog
        %tstress
        %tstrain
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
            
            obj.stress = zeros(dim.nstre,nelem,G.ngaus);
            obj.strain = zeros(dim.nstre,nelem,G.ngaus);
            
             obj.stress_homog = zeros(dim.nstre,1);
            
            for igaus=1:G.ngaus
                vol_dom = sum(G.dvolu(:,igaus));

                obj.strain(1:dim.nstre,:,igaus) = vstrain'*ones(1,nelem) + obj.strain_fluct(1:dim.nstre,:,igaus);
                for istre=1:dim.nstre
                    for jstre=1:dim.nstre
                        obj.stress(istre,:,igaus) = squeeze(obj.stress(istre,:,igaus)) + 1/vol_dom*squeeze(squeeze(Cmat(istre,jstre,:)))'.* obj.strain(jstre,:,igaus);
                    end
                end
                % contribucion a la C homogeneizada
                for istre=1:dim.nstre
                    obj.stress_homog(istre) = obj.stress_homog(istre) +  1/vol_dom *(obj.stress(istre,:,igaus))*G.dvolu(:,igaus);
                end                
            end
        end
    end
end

