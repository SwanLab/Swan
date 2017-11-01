classdef Element<handle
    %Element Summary of this class goes here
    %   Detailed explanation goes here TEST
    
    properties (GetAccess = ?Physical_Problem, SetAccess = protected)
        RHS
        LHS
    end
    
    properties (GetAccess = {?Element_Elastic,?Element_Thermal,?PhysicalVariables}, SetAccess = {?Physical_Problem,?Element})
        B
    end
    
    methods     
    end
    
end
