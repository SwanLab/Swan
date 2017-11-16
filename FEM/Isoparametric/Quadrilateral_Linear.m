classdef Quadrilateral_Linear < Isoparametric
    properties
    end
    
    methods
        function obj = Quadrilateral_Linear
            obj = obj@Isoparametric;
            obj.type = 'QUAD';
            obj.ndime = 2;
            obj.nnode = 4;
            obj.ngaus = 1;
            
            % Compute WEIGP and POSGP
            obj.posgp = [0,0];
            obj.weigp = 1;
            
            for igaus = 1:obj.ngaus
                s = obj.posgp(1,igaus);
                t = obj.posgp(2,igaus);
                st = s*t;
                
                % Shape Functions
                obj.shape(1,igaus) = (1.-t-s+st)*0.25;
                obj.shape(2,igaus) = (1.-t+s-st)*0.25;
                obj.shape(3,igaus) = (1.+t+s+st)*0.25;
                obj.shape(4,igaus) = (1.+t-s-st)*0.25;
                
                % SF Derivatives
                obj.deriv(1,1,igaus) = (-1.+t)*0.25;
                obj.deriv(1,2,igaus) = (+1.-t)*0.25;
                obj.deriv(1,3,igaus) = (+1.+t)*0.25;
                obj.deriv(1,4,igaus) = (-1.-t)*0.25;
                obj.deriv(2,1,igaus) = (-1.+s)*0.25;
                obj.deriv(2,2,igaus) = (-1.-s)*0.25;
                obj.deriv(2,3,igaus) = (+1.+s)*0.25;
                obj.deriv(2,4,igaus) = (+1.-s)*0.25;
            end
        end
    end
    
end