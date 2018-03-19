classdef Material_Hyperelastic < Material_Elastic
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Element_Hyperelastic, ?Material_Hyperelastic_2D, ?Material_Hyperelastic_3D, ?PhysicalVars_Elastic, ?Physical_Problem}, SetAccess = ?Physical_Problem)
%         connec
%         cartd0
%         nnode
%         coord
    end
    
    methods
        function obj = Material_Hyperelastic(nelem)
            obj@Material_Elastic(nelem);
%             obj.connec= connec;
%             obj.cartd0 = cartd;
%             obj.nnode = nnode;
%             obj.coord = coord;
        end
    end
end