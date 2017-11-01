classdef Material_Elastic_ISO < Material_Elastic
    %Material_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Material_Elastic_ISO_2D, ?Material_Elastic_ISO_3D, ?PhysicalVars_Elastic}, SetAccess = protected)
        kappa
        mu
    end
    
    % !! Property attributes will have to be changed when OPT_problem is implemented !!
    methods (Access = protected)
        function obj = Material_Elastic_ISO(nelem)
            obj.nelem = nelem;
            obj.kappa = .9107;
            obj.mu = .3446;
        end
    end
    
    methods (Access = protected)
        function obj = setProps_parent(props)
            obj.kappa = props.kappa;
            obj.mu = props.mu;
        end
    end
end

