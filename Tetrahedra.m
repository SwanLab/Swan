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
            obj.ngaus = 1;          % tetrahedra
            obj.weigp = 1/6;
            obj.posgp = [1/4;1/4;1/4];
            
            % s : xi coordinate
            % t : eta coordinate
            % u : zeta coordinate (for 3D)
            s = obj.posgp(1,obj.ngaus);
            t = obj.posgp(2,obj.ngaus);
            u = obj.posgp(3,obj.ngaus);
            obj.deriv = zeros(obj.ndime,obj.nnode);
            obj.shape = zeros(1,obj.nnode);
            
            % Shape Functions
            obj.shape(1)=(1-t-s-u);                    
            obj.shape(2)=s;           
            obj.shape(3)=t;           
            obj.shape(4)=u;           
       
            
            % Derivatives
            obj.deriv=[-1 1 0 0
                             -1 0 1 0
                             -1 0 0 1];
        end
    end
    
end
