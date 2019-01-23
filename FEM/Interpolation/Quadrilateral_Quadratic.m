classdef Quadrilateral_Quadratic < Interpolation
    
    methods (Access = public)
        function obj = Quadrilateral_Quadratic
            obj = obj@Interpolation;
            obj.type = 'QUAD';
            obj.ndime = 2;
            obj.nnode = 8;
            obj.pos_nodes = [-1 -1; 1 -1; 1 1; -1 1];
            
            shape = @(s,t) {(1.-t-s+s*t)*0.25;
                            (1.-t+s-s*t)*0.25;
                            (1.+t+s+s*t)*0.25;
                            (1.+t-s-s*t)*0.25};
            obj.shape = shape;
            
            deriv = @(s,t) {(-1.+t)*0.25 (+1.-t)*0.25 (+1.+t)*0.25 (-1.-t)*0.25;
                            (-1.+s)*0.25 (-1.-s)*0.25 (+1.+s)*0.25 (+1.-s)*0.25};
            obj.deriv = deriv;
            % Compute WEIGP and POSGP
%             a =  0.577350269189626;
%             obj.posgp(1,:) = [-a,-a];
%             obj.posgp(2,:) = [+a,-a];
%             obj.posgp(3,:) = [-a,+a];
%             obj.posgp(4,:) = [+a,+a];
%             obj.weigp = 0.25*ones(1,ngaus);
%             
%             for igaus = 1:obj.ngaus
%                 s = obj.posgp(1,igaus);
%                 t = obj.posgp(2,igaus);
%                 st = s*t;
%                 
%                 % Shape Functions
%                 obj.shape(1,igaus) = (1.-t-s+st)*0.25;
%                 obj.shape(2,igaus) = (1.-t+s-st)*0.25;
%                 obj.shape(3,igaus) = (1.+t+s+st)*0.25;
%                 obj.shape(4,igaus) = (1.+t-s-st)*0.25;
%                 
%                 % SF Derivatives
%                 obj.deriv(1,1,igaus) = (-1.+t)*0.25;
%                 obj.deriv(1,2,igaus) = (+1.-t)*0.25;
%                 obj.deriv(1,3,igaus) = (+1.+t)*0.25;
%                 obj.deriv(1,4,igaus) = (-1.-t)*0.25;
%                 obj.deriv(2,1,igaus) = (-1.+s)*0.25;
%                 obj.deriv(2,2,igaus) = (-1.-s)*0.25;
%                 obj.deriv(2,3,igaus) = (+1.+s)*0.25;
%                 obj.deriv(2,4,igaus) = (+1.-s)*0.25;
%             end
        end
    end
    
end