classdef Quadrature_Triangle < Quadrature

    methods
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case 'CONSTANT'
                    obj.ngaus = 1;
                    obj.weigp = 1/2;    
                    obj.posgp = [1/3;1/3];
                    
                case 'LINEAR'
                    obj.ngaus = 1;        
                    obj.weigp = [1/2];
                    obj.posgp = [1/3;1/3];
                    
                case 'QUADRATIC'
                    %obj.ngaus = 3;         
                    %obj.weigp = [1/6;1/6;1/6];
                    %obj.posgp = [0,0.5;0.5,0;0.5,0.5]';
                    obj.ngaus = 3;         
                    obj.weigp = [1/6;1/6;1/6];
                    obj.posgp = [2/3,1/6;1/6,1/6;1/6,2/3]';                    
                    
                case 'QUADRATIC2'
                    obj.ngaus = 2;         
                    obj.weigp = [1/6;1/6;1/6];
                    obj.posgp = [1/6,1/6;2/3,1/6;1/6,2/3]';                    
                    
                case 'QUADRATICMASS'
                    obj.ngaus = 3;          
                    obj.weigp = [1/6;1/6;1/6];
                    obj.posgp = [2/3,1/6,1/6;1/6,2/3,1/6];
                otherwise
                    disp('Quadrature not implemented for triangle elements')
            end
        end
    end
end