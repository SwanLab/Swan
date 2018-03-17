classdef Tetrahedra<Interpolation
    
    properties
    end
    
    methods
        % Constructor
        function obj = Tetrahedra(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'TETRAHEDRA';
            obj.order= 'LINEAR';
            obj.ndime = 3;
            obj.nnode = 4;
        end
        function computeShapeDeriv(obj,posgp)
            for igaus=1:size(posgp,2)
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                u = posgp(3,igaus);
                obj.shape(:,igaus) =[(1.-t-s-u);
                    s;
                    t;
                    u];
                
                % Derivatives
                obj.deriv(:,:,igaus)=[-1 1 0 0
                    -1 0 1 0
                    -1 0 0 1];
            end
        end
    end
    
end
