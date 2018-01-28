classdef Element<handle
    %Element Summary of this class goes here
    %   Detailed explanation goes here TEST
    
    %% !! NEEDS REVISION !! -> should B be a class?? Or just be contained in element ??
    
    properties (GetAccess = {?Physical_Problem, ?Element_Elastic_Micro}, SetAccess = protected)
        RHS
        LHS
    end
    
    properties (GetAccess = {?Element_Elastic,?Element_Thermal,?PhysicalVariables,?Element_Elastic_Micro}, SetAccess = {?Physical_Problem,?Element, ?Element_Elastic_Micro})
        B
    end
    
    methods (Access = ?Physical_Problem, Static)
        function element = create(ptype,pdim)
            switch ptype
                case 'ELASTIC'
                    switch pdim
                        case '2D'
                            element = Element_Elastic;
                            element.B = B2;
                        case '3D'
                            element = Element_Elastic;
                            element.B = B3;
                    end
                case 'THERMAL'
                    element = Element_Thermal;
                    element.B = B_thermal;
                otherwise
                    error('Invalid ptype.')
            end
        end
    end
    methods (Access = {?Physical_Problem, ?Element})
        function obj = computeRHS(obj,nunkn,nelem,nnode,bc,idx)
            RHSSuperficial  = obj.computeSuperficialRHS(nunkn,nelem,nnode,bc,idx);
            RHSVolumetric  = obj.computeVolumetricRHS(nunkn,nelem,nnode,bc,idx);
            obj.RHS = RHSSuperficial + RHSVolumetric;
        end
    end
    
    methods (Abstract, Access = protected)
        r = computeSuperficialRHS(obj)
        r = computeVolumetricRHS(obj)
    end
end

