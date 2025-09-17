classdef RigidBodyModesComputer < handle

    properties (Access = private)
        mesh
    end

    methods (Access = public)
        function obj = RigidBodyModesComputer(m)
            obj.mesh = m;
        end

        function R = compute(obj,refPoint)
            ndim = obj.mesh.ndim;

            % SWITCH
            
   % rigidBodyModes3D computes the 6 rigid body modes in 3D
            %
            % INPUTS:
            %   coords    : N x 3 matrix of nodal coordinates [x, y, z]
            %   refPoint  : 1 x 3 vector [x0, y0, z0] - reference point for rotation
            %
            % OUTPUT:
            %   R         : 3N x 6 matrix where each column is a rigid body mode
            %               (3 displacements per node)

            coords = obj.mesh.coord;
            % Number of nodes
            N = size(coords, 1);

            % Initialize rigid body modes matrix
            R = zeros(3*N, 6);

            % Displacement indices
            ix = 1:3:3*N;
            iy = 2:3:3*N;
            iz = 3:3:3*N;

            % Displacements for translation modes
            R(ix, 1) = 1; % Translation in x
            R(iy, 2) = 1; % Translation in y
            R(iz, 3) = 1; % Translation in z

            % Compute position vectors relative to reference point
            relCoords = coords - refPoint;

            x = relCoords(:, 1);
            y = relCoords(:, 2);
            z = relCoords(:, 3);

            % Displacements for rotation about x-axis
            R(iy, 4) = -z;
            R(iz, 4) =  y;

            % Displacements for rotation about y-axis
            R(ix, 5) =  z;
            R(iz, 5) = -x;

            % Displacements for rotation about z-axis
            R(ix, 6) = -y;
            R(iy, 6) =  x;
        end
    end
end