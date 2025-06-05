classdef LHSIntegratorNavierStokes < handle

    properties (Access = private)
        mesh
        velocityFun
        pressureFun
        velocityField
        material
        LHSS
        dt
    end

    methods (Access = public)

        function obj = LHSIntegratorNavierStokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            C   = obj.computeConvectiveMatrix();
            LHS = obj.addCMToLHS(C);
        end

        function C = computeConvectiveTerm(obj)
            C   = obj.computeConvectiveMatrix();
        end

        function [M,D,KU,KP] = computeLHSMatrix(obj)
            M  = obj.computeMassMatrix();
            D  = obj.computeWeakDivergenceMatrix();
            KU = obj.computeVelocityLaplacian();       
            KP = obj.computePressureMatrix();
        end

    end

    methods (Access = private)
    
        function init(obj, cParams)
            obj.mesh          = cParams.mesh;
            obj.material      = cParams.material;
            obj.velocityFun   = cParams.velocityFun;
            obj.pressureFun   = cParams.velocityFun;
            obj.velocityField = cParams.velocityField;
            obj.LHSS          = cParams.LHSS;
            obj.dt            = cParams.dt;
        end

        function c = computeConvectiveMatrix(obj)
            s.type  = 'Convective';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.material = obj.material;
            C = LHSIntegrator.create(s);
            c = C.compute(obj.velocityField);
        end

        function D = computeWeakDivergenceMatrix(obj)
            s.type = 'WeakDivergence';
            s.mesh = obj.mesh;
            s.trial = obj.pressureFun;
            s.test  = obj.velocityFun;
            LHS = LHSIntegrator.create(s);
            D = LHS.compute();
        end

        function lhs = computeVelocityLaplacian(obj)
            s.type  = 'Laplacian';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.material = obj.material;
            LHS = LHSIntegrator.create(s);
            lhs = LHS.compute();
            lhs = obj.symGradient(lhs);
        end
        
        function A = symGradient(obj, B)
            A = 1/2 * (B+B');
        end

        function M = computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.quadratureOrder = 3;
            LHS = LHSIntegrator.create(s);
            m = LHS.compute();

            dtime = obj.dt;
            M = m/dtime;
        end

        function KP = computePressureMatrix(obj)
            s.type  = 'StiffnessMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.pressureFun;
            s.trial = obj.pressureFun;
            LHS = LHSIntegrator.create(s);
            KP  = LHS.compute();

        end

        function LHS = addCMToLHS(obj,C)
            LHS  = obj.LHSS;
            nVel = size(C, 1);
            LHS(1:nVel, 1:nVel) = LHS(1:nVel, 1:nVel) + C;
        end

    end

end