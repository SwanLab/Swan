classdef Tetrahedra<Isoparametric
    
    properties
    end
    
    methods
        % Constructor
        function obj = Tetrahedra()
            obj = obj@Isoparametric();
            obj.type = 'TETRAHEDRA';
            obj.ndime = 3;          % 1D/2D/3D
            obj.nnode = 4;
%             obj.ngaus = 1;          % tetrahedra
%             obj.weigp = 1/6;
%             obj.posgp = [1/4;1/4;1/4];           
            % s : xi coordinate
            % t : eta coordinate
            % u : zeta coordinate (for 3D)

            shape =@(s,t,u) {(1.-t-s-u);
                                s;
                                t;
                                u};
            obj.shape = shape;    
            % Derivatives
            deriv=@(x,y,z){[-1 1 0 0
                             -1 0 1 0
                             -1 0 0 1]};
            obj.deriv=deriv;
        end
    end
    
end
