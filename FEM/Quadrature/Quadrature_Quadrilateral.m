classdef Quadrature_Quadrilateral<Quadrature
    properties
    end
    methods
        function computeQuadrature(obj,order)
            switch order
                case 'LINEAR'
                    obj.ngaus = 4;                    
                    % Compute WEIGP and POSGP
                    a =  0.577350269189626;
                    obj.posgp(:,1) = [-a,-a];
                    obj.posgp(:,2) = [+a,-a];
                    obj.posgp(:,3) = [-a,+a];
                    obj.posgp(:,4) = [+a,+a];
                    
                    obj.weigp =  [1,1,1,1];%1*ones(1,obj.ngaus);
                case {'QUADRATIC','QUADRATICMASS'} %SERENDIPITY, QUADRILATERAL QUADRATIC NOT IMPLEMENTED
                    obj.ngaus = 9;                    
                    % Compute WEIGP and POSGP
                    a =  0.77459667;
                    obj.posgp(:,1) = [ 0,+a];
                    obj.posgp(:,2) = [ 0, 0];
                    obj.posgp(:,3) = [+a,+a];
                    obj.posgp(:,4) = [-a,-a];
                    obj.posgp(:,5) = [-a, 0];
                    obj.posgp(:,6) = [+a, 0];
                    obj.posgp(:,7) = [+a,-a];
                    obj.posgp(:,8) = [-a,+a];
                    obj.posgp(:,9) = [ 0,-a];
                    
                    obj.weigp =ones(1,obj.ngaus);
                otherwise
                    error('Invalid interpolation order for element Quadrilateral.');
            end
        end
    end
end