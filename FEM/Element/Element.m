classdef Element<handle
    %Element Summary of this class goes here
    %   Detailed explanation goes here TEST
    
    %% !! NEEDS REVISION !! -> should B be a class?? Or just be contained in element ??
    
    properties (GetAccess = ?Physical_Problem, SetAccess = protected)
        RHS
        LHS
    end
    
    properties (GetAccess = {?Element_Elastic,?Element_Thermal,?PhysicalVariables}, SetAccess = {?Physical_Problem,?Element})
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
    methods (Access = ?Physical_Problem)
        function obj = computeRHS(obj,nunkn,nelem,nnode,bc,idx)
            RHSPuntual = obj.computePuntualRHS(nunkn,nelem,nnode,bc,idx);
            RHSSuperficial  = obj.computeSuperficialRHS(nunkn,nelem,nnode,bc,idx);
            RHSVolumetric  = obj.computeVolumetricRHS(nunkn,nelem,nnode,bc,idx);
            obj.RHS = RHSSuperficial + RHSVolumetric + RHSPuntual;
        end
    end
    
    methods (Abstract, Access = protected)
        r = computePuntualRHS(obj)
        r = computeSuperficialRHS(obj)
        r = computeVolumetricRHS(obj)
    end
end

