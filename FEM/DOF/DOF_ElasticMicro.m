classdef DOF_ElasticMicro < DOF_Elastic
    %DOF_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = DOF_ElasticMicro(problemID,geometry,mesh)
            obj@DOF_Elastic(problemID,geometry,mesh);
        end
        
        function constrained = compute_constrained_dof(obj,ifield)
            constrained = [obj.periodic_constrained;obj.dirichlet{ifield}];
        end
    end
end

