classdef Quadrature < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        posgp
        weigp
        ngaus
        order
    end
    methods
        function computeQuadrature(obj,order)
            obj.posgp=[];
            obj.weigp=[];
            obj.order=order;
        end
    end
    methods (Static)        
        function quadrature=set(type)
            switch type
                case 'LINE'
                    quadrature = Quadrature_Line;
                case 'TRIANGLE'
                    quadrature = Quadrature_Triangle;                    
                case 'QUAD'                    
                    quadrature = Quadrature_Quadrilateral;                    
                case 'TETRAHEDRA'
                    quadrature = Quadrature_Tetrahedra;
                case 'HEXAHEDRA'
                    quadrature = Quadrature_Hexahedra;
                otherwise
                    error('Invalid quadrature type.')
            end
        end
    end
end

