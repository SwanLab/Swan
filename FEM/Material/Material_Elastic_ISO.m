classdef Material_Elastic_ISO < Material_Elastic
    %Material_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    

    properties (GetAccess = public, SetAccess = protected)
        kappa
        mu
        lambda
    end
    
    % !! Property attributes will have to be changed when OPT_problem is implemented !!
    methods (Access = protected)
        function obj = Material_Elastic_ISO(nelem)
            obj@Material_Elastic(nelem);

            obj.kappa  = 1;.9107;
            obj.mu     = 0.5;.3446;            
            obj.lambda = obj.kappa-obj.mu;
            obj = obj.computeC;
        end
    end
    
    methods %(Access = ?Physical_Problem)
        function obj = setProps(obj,props)
            obj.kappa = props.kappa;
            obj.mu = props.mu;
            obj.lambda = obj.kappa-obj.mu;
            obj.computeC;
        end
    end
    
    methods (Access = protected)
        function obj = computeC(obj)
        end
    end
end

