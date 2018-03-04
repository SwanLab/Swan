classdef Triangle_Quadratic<Isoparametric
    %Triangle_Quadratic Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Constructor
        function obj = Triangle_Quadratic
            obj = obj@Isoparametric;
            obj.type = 'TRIANGLE_QUADRATIC';
            obj.ndime = 2;          % 1D/2D/3D
            obj.nnode = 6;
            obj.ngaus = 3;          % Linear triangle
            obj.weigp = [1/3;1/3;1/3];
            obj.posgp = [0,0.5;0.5,0;0.5,0.5]';
            obj.pos_nodes = [0,0 ; 1 0; 0,1 ; 0.5,0 ; 0.5,0.5 ; 0,0.5];
%             [0,0;0.5 0;1,0; 0.5 0.5; 0,1; 0 0.5]
%             obj.posgp = [1/6,1/6;2/3,1/6;1/6,2/3]';
            
            % s : xi coordinate
            % t : eta coordinate
            % u : zeta coordinate (for 3D)
            s = obj.posgp(1,1:obj.ngaus);
            t = obj.posgp(2,1:obj.ngaus);
            obj.deriv = zeros(obj.ndime,obj.nnode);
            
            % !! Originally: obj.shape = zeros(1,obj.nnode);
            shape = @(s,t) {(1.0-s-t)*(1.0-2*s-2*t),s*(2*s-1.0),t*(2*t-1.0),4*s*(1.0-s-t),4*s*t,4*t*(1.0-s-t)};
            
            obj.shape = shape;
            % Shape Functions
%             obj.shape(1,:) = (1.0-s-t).*(1.0-2*s-2*t);
%             obj.shape(2,:) = s.*(2.*s-1.0);
%             obj.shape(3,:) = t.*(2.*t-1.0);
%             obj.shape(4,:) = 4.*s.*(1.0-s-t);
%             obj.shape(5,:) = 4.*s.*t;
%             obj.shape(6,:) = 4.*t.*(1.0-s-t);
            
            % Derivatives
            %% !! NEEDS REVISION (VERIFY DERIVATIVES) !! (http://www.ce.memphis.edu/7111/notes/class_notes/chapter_03d_slides.pdf)
            % w.r.t. xi
            
            deriv = @(s,t) {4*(s+t)-3    4*s-1  0.0     4*(1.0-t)-8*s   4*t     -4*t;
                            4*(s+t)-3.0  0.0    4*t-1.0     -4*s        4*s     4*(1.0-s)-8*t};
            obj.deriv=deriv;
%             obj.deriv(1,1,1:obj.ngaus) = 4.*(s+t)-3;
%             obj.deriv(1,2,1:obj.ngaus) = 4.*s-1;
%             obj.deriv(1,3,1:obj.ngaus) = 0.0;
%             obj.deriv(1,4,1:obj.ngaus) = 4.*(1.0-t)-8.*s;
%             obj.deriv(1,5,1:obj.ngaus) = 4.*t;
%             obj.deriv(1,6,1:obj.ngaus) = -4.*t;
%             % w.r.t. eta
%             obj.deriv(2,1,1:obj.ngaus) = 4.*(s+t)-3.0;
%             obj.deriv(2,2,1:obj.ngaus) = 0.0;
%             obj.deriv(2,3,1:obj.ngaus) = 4.*t-1.0;
%             obj.deriv(2,4,1:obj.ngaus) = -4.*s;
%             obj.deriv(2,5,1:obj.ngaus) = 4.*s;
%             obj.deriv(2,6,1:obj.ngaus) = 4.*(1.0-s)-8.*t;            
        end   
    end
end