classdef PhysicalVars_Elastic_3D < PhysicalVars_Elastic
    %PhysicalVars_Elastic_3D Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! STUDY TO INTRODUCE JUST C !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function obj = computeVars(obj,d_u,dim,nnode,nelem,ngaus,idx,element,material)
            nstre = 6;
            obj.d_u = d_u;
            strain = obj.computeStrain(d_u,dim,nnode,nelem,ngaus,idx,element);
            obj.strain = strain;
            obj.stress = obj.computeStress(strain,material.C,ngaus,nstre);
        end
    end
    
end

