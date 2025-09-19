classdef RigidBodyModesComputer < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        function obj = RigidBodyModesComputer(m)
            obj.mesh = m;
        end

        function R = compute(obj,refPoint)
            coords    = obj.mesh.coord;
            relCoords = coords - refPoint;
            switch obj.mesh.ndim
                case 3
                    x = relCoords(:, 1);
                    y = relCoords(:, 2);
                    z = relCoords(:, 3);
                    R = obj.compute3D(x,y,z);
            end
        end
    end

    methods (Static, Access = private)
        function R = compute3D(x,y,z)
            N  = length(x);
            R  = zeros(3*N, 6);
            ix = 1:3:3*N;
            iy = 2:3:3*N;
            iz = 3:3:3*N;

            R(ix, 1) = 1;
            R(iy, 2) = 1;
            R(iz, 3) = 1;
            
            R(iy, 4) = -z;
            R(iz, 4) =  y;
            R(ix, 5) =  z;
            R(iz, 5) = -x;
            R(ix, 6) = -y;
            R(iy, 6) =  x;
        end
    end
end