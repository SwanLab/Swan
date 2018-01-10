
classdef PhysicalVars_Elastic_2D_Micro < PhysicalVars_Elastic_2D
    %PhysicalVars_Elastic_2D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties (Access = {?Physical_Problem_Micro, ?PhysicalVars_Elastic_2D})
        stress_homog
    end
    
    methods (Access = {?Physical_Problem, ?PhysicalVars_Elastic_2D, ?PhysicalVars_Elastic_2D_Micro})
        function obj = PhysicalVars_Elastic_2D_Micro(ndof)
            obj.d_u = zeros(ndof,1);      
        end
        
        function obj = computeVars(obj,d_u,dim,G,nelem,idx,element,material,nstre)
            obj = computeVars@PhysicalVars_Elastic_2D(obj,d_u,dim,G,nelem,idx,element,material,nstre);
            stress_h = zeros(dim.nstre,1);
            for ielem=1:nelem
                for igaus=1:G.ngaus
                    for istres=1:nstre
                        stress_h(istres) = stress_h(istres) + obj.stress(igaus,istres,ielem)*G.dvolu(ielem,igaus);
                    end
                end
            end
            obj.stress_homog = stress_h / sum(G.dvolu);
        end
    end
end

