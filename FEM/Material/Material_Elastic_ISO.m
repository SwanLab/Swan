classdef Material_Elastic_ISO < Material_Elastic
    %Material_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Material_Elastic_ISO_2D, ?Material_Elastic_ISO_3D, ?PhysicalVars_Elastic}, SetAccess = protected)

    end
    
    % !! Property attributes will have to be changed when OPT_problem is implemented !!
    methods (Access = protected)
        function obj = Material_Elastic_ISO(nelem)
            obj@Material_Elastic(nelem);
            obj = obj.computeC;
        end
    end
    
    methods (Access = ?Physical_Problem)
        function obj = setProps(obj,props)
            obj.kappa = props.kappa;
            obj.mu = props.mu;
            obj = obj.computeC;
        end
    end
    
    methods (Access = protected)
        function obj = computeC(obj)
        end
    end
end

