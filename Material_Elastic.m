classdef Material_Elastic < Material
    %Material_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! This has to be revised. !!
    % !! The best structure will depend on how this is wanted to be initialized !!
    
    properties (GetAccess = {?Element_Elastic,?PhysicalVars_Elastic}, SetAccess = protected)
        C
    end
    
    methods
    end    
end

