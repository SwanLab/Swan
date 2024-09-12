classdef ComputeInitialEnergy < handle

    properties (Access = public)
        e0
    end

    properties (Access = private)
        mesh
        bc
        C
    end

    methods (Access = public)

        function obj = ComputeInitialEnergy(cParams)
            obj.init(cParams)
            obj.computeEnergy();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.C = cParams.C;
            obj.bc = cParams.bc;
            obj.mesh = cParams.mesh;
        end

        function computeEnergy(obj)
            [U,F] = obj.solveFEM();

            obj.e0 = 0.5 * dot(F,U);
        end

        function [U,F] = solveFEM(obj)
            s.mesh               = obj.mesh;
            s.scale              = 'MACRO';
            s.material           = obj.C;
            s.dim                = '2D';
            s.boundaryConditions = obj.bc;
            s.solverType         = 'REDUCED';
            s.solverMode         = 'DISP';

            fem         = ElasticProblem(s);
            fem.solve();
            displ    = fem.uFun.fValues; 
            %obj.displFun = fem.uFun;
            U          = reshape(displ, [size(displ,1)*2, 1]);
            force      = fem.forces; 
            Fx         = force(1:2:end); % Fx values are at odd indices
            Fy         = force(2:2:end); % Fy values are at even indices
            forcesVect = [Fx Fy];
            F          = reshape(forcesVect, [size(displ,1)*2, 1]);
        end
    end
end