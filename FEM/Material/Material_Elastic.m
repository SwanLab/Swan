classdef Material_Elastic < Material
    %Material_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! This has to be revised. !!
    % !! The best structure depends on how this is wanted to be initialized !!
    
    properties (GetAccess = {?Element,?Material_Elastic_ISO,?Material_Hyperelastic,?Physical_Problem}, SetAccess = public) 
        C
        kappa
        mu
        lambda
    end
    
    methods (Access = protected)
        function obj = Material_Elastic(nelem)
            obj@Material(nelem);
            obj.kappa  = .9107;
            obj.mu     = .3446;
            obj.lambda = obj.kappa-obj.mu;
        end
    end
end
