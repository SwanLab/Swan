classdef LHSintegrator_Stokes < handle %LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
    end

    properties (Access = private)
        dt
        mesh
        velocityField
        pressureField
        material
        D
    end

    methods (Access = public)

        function obj = LHSintegrator_Stokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            velLHS = obj.computeVelocityLHS();
            D      = obj.computeDmatrix();
            prsLHS = obj.computePressureLHS(D);
            LHS = [velLHS, D; D',prsLHS];
        end

    end

    methods (Access = private)
    
        function init(obj, cParams)
            obj.dt            = cParams.dt;
            obj.mesh          = cParams.mesh;
            obj.material      = cParams.material;
            obj.pressureField = cParams.pressureField;
            obj.velocityField = cParams.velocityField;
        end

        function LHS = computeVelocityLHS(obj)
            K = obj.computeVelocityLaplacian();
            M = obj.computeMassMatrix();
            lhs = K + M;
            LHS = obj.symGradient(lhs);
        end

        function D = computeDmatrix(obj)
            s.type = 'StokesD';
            s.mesh = obj.mesh;
            s.pressure = obj.pressureField;
            s.velocity = obj.velocityField;
            LHS = LHSintegrator.create(s);
            D = LHS.compute();
        end

        function BB = computePressureLHS(obj,D)
            sz = size(D, 2);
            BB = sparse(sz,sz);
        end

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end

        function lhs = computeVelocityLaplacian(obj)
            s.type  = 'Laplacian';
            s.mesh  = obj.mesh;
            s.field = obj.velocityField;
            s.material = obj.material;
            LHS = LHSintegrator.create(s);
            lhs = LHS.compute();
            lhs = obj.symGradient(lhs);
        end

        function M = computeMassMatrix(obj)
            vel = obj.velocityField;
            dtime = obj.dt;
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = vel;
            LHS = LHSintegrator.create(s);
            m = LHS.compute();
            M = m/dtime;
            obj.M = M;
        end

    end

end