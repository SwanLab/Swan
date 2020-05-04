classdef Quadrature_Line < Quadrature

    methods
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case 'CONSTANT'
                    obj.ngaus = 1;
                    obj.weigp = 2;    
                    obj.posgp = [0];
                    
                case 'LINEAR'
                    obj.ngaus = 2;        
                    obj.weigp = [1,1];
                    obj.posgp = [-1/sqrt(3),1/sqrt(3)];
                    
                case 'QUADRATIC'
                    obj.ngaus = 3;         
                    obj.weigp = [5/9,8/9,5/9];
                    obj.posgp = [-sqrt(3/5),0,sqrt(3/5)];
                    
                otherwise
                    disp('Quadrature not implemented for triangle elements')
            end
        end
    end
end