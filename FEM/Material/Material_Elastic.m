classdef Material_Elastic < Material
    %Material_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! This has to be revised. !!
    % !! The best structure depends on how this is wanted to be initialized !!
    
    properties (GetAccess = {?Element,?Material_Elastic_ISO,?Material_Hyperelastic_2D,?PhysicalVars_Elastic}, SetAccess = protected) 
        C
    end
    
    methods (Access = protected)
        function obj = Material_Elastic(nelem)
            obj@Material(nelem);
        end
    end
end
