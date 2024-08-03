classdef LHSintegrator_Stokes < handle %LHSintegrator

    properties (GetAccess = public, SetAccess = private)
        M
    end

    properties (Access = private)
        dt
        mesh
        velocityFun
        pressureFun
        material
    end

    methods (Access = public)

        function obj = LHSintegrator_Stokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            velLHS = obj.computeVelocityLHS();
            D      = obj.computeWeakDivergenceMatrix();
            prsLHS = obj.computePressureLHS();
            LHS = [velLHS, D; D',prsLHS];
        end

    end

    methods (Access = private)
    
        function init(obj, cParams)
            obj.dt          = cParams.dt;
            obj.mesh        = cParams.mesh;
            obj.material    = cParams.material;
            obj.pressureFun = cParams.pressureFun;
            obj.velocityFun = cParams.velocityFun;
        end

        function LHS = computeVelocityLHS(obj)
            K = obj.computeVelocityLaplacian();
            M = obj.computeMassMatrix();
            lhs = K + M;
            LHS = obj.symGradient(lhs);
        end

        function D = computeWeakDivergenceMatrix(obj)
            s.type = 'WeakDivergence';
            s.mesh = obj.mesh;
            s.trial = obj.velocityFun;
            s.test  = obj.pressureFun;
            LHS = LHSintegrator.create(s);
            D = LHS.compute();
        end

        function BB = computePressureLHS(obj)
            ndofs = obj.pressureFun.nDofs;
            BB = sparse(ndofs,ndofs);
        end

        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end

        function lhs = computeVelocityLaplacian(obj)
            s.type  = 'Laplacian';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.material = obj.material;
            LHS = LHSintegrator.create(s);
            lhs = LHS.compute();
            lhs = obj.symGradient(lhs);
        end

        function M = computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.quadratureOrder = 3;
            LHS = LHSintegrator.create(s);
            m = LHS.compute();

            dtime = obj.dt;
            M = m/dtime;
            obj.M = M;
        end

    end

end