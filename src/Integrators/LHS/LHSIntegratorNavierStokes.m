classdef LHSIntegratorNavierStokes < handle

    properties (Access = private)
        mesh
        velocityFun
        pressureFun
        velocityField
        material
        LHSS
        dt
        residual
        nu
    end

    methods (Access = public)

        function obj = LHSIntegratorNavierStokes(cParams)
            obj.init(cParams);
        end

        function LHS = compute(obj)
            [C, ~] = obj.computeConvectiveMatrix();
            LHS = obj.addCMToLHS(C);
        end

        function LHS = computeStabilization(obj)
            [C , a]   = obj.computeConvectiveMatrix();
            Kvva2   = obj.computeVVA2Mat();
            Kavp    = obj.computeAVPMat();
            Kpva    = Kavp';
            Kpp     = obj.computePPMat();
            Kpdv    = obj.computePDivergenceVMatrix();
            Kdvdv   = obj.computeVDVDMatrix();
            tauS1 = obj.computeTauS1(C, Kvva2);
            tauS3 = obj.computeTauS3(tauS1, a);
            tauS   = obj.computeTau(tauS1,tauS3);
            tauP1 = obj.computeTauS1(Kpdv, Kpva);
            tauP3 = obj.computeTauP3(tauP1, a, tauS1);
            tauP   = obj.computeTau(tauP1,tauP3);
            tauL   = obj.computeTauS1(C,Kdvdv);
            LHS   = obj.addTermsToLHS(C, tauS, Kvva2, Kavp, tauP, Kpva, Kpp, tauL, Kdvdv);
        end



        function [C , a] = computeConvectiveTermAndVelocityF(obj)
            [C , a]  = obj.computeConvectiveMatrix();
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
            obj.pressureFun   = cParams.pressureFun;
            obj.velocityField = cParams.velocityField;
            obj.LHSS          = cParams.LHSS;
            obj.dt            = cParams.dt;
            obj.residual      = cParams.residual;
            obj.nu            = cParams.material.nuValue;
        end

        function [c,a] = computeConvectiveMatrix(obj)
            s.type  = 'Convective';
            s.mesh  = obj.mesh;
            s.test  = obj.velocityFun;
            s.trial = obj.velocityFun;
            s.material = obj.material;
            C = LHSIntegrator.create(s);
            [c,a] = C.compute(obj.velocityField);
        end

        function D = computeWeakDivergenceMatrix(obj)
            s.type = 'WeakDivergence';
            s.mesh = obj.mesh;
            s.trial = obj.pressureFun;
            s.test  = obj.velocityFun;
            LHS = LHSIntegrator.create(s);
            D = LHS.compute();
        end

        function D = computePDivergenceVMatrix(obj)
            s.type = 'WeakDivergence';
            s.mesh = obj.mesh;
            s.trial = obj.pressureFun;
            s.test  = obj.velocityFun;
            LHS     = LHSIntegrator.create(s);
            D       = LHS.compute();
            D       = - D';
        end

        function DD = computeVDVDMatrix(obj)
            s.type = 'DobleDivergenceMatrix';
            s.mesh = obj.mesh;
            s.trial = obj.velocityFun;
            s.test  = obj.velocityFun;
            LHS     = LHSIntegrator.create(s);
            DD       = - LHS.compute();
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

        function Kvs = computeVVA2Mat(obj)
            s.type = 'StiffnessMatrixWithFunction';
            s.mesh       = obj.mesh;
            s.trial      = obj.velocityFun;
            s.test       = obj.velocityFun;
            s.function   = obj.velocityField;
            s.function.setFValues(s.function.fValues .^2);
            LHS = LHSIntegrator.create(s);
            Kvs = LHS.compute();
        end

        function K = computeAVPMat(obj)
            s.type       = 'StiffnessMatrixWithFunction';
            s.mesh       = obj.mesh;
            s.trial      = obj.pressureFun;
            s.test       = obj.velocityFun;
            s.function   = obj.velocityField;
            LHS = LHSIntegrator.create(s);
            K = LHS.compute();
        end

        function K = computePPMat(obj)
            s.type       = 'StiffnessMatrix';
            s.mesh       = obj.mesh;
            s.trial      = obj.pressureFun;
            s.test       = obj.pressureFun;
            LHS = LHSIntegrator.create(s);
            K = LHS.compute();
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

        function LHS = addTermsToLHS(obj, C, tauS, Kvva2, Kavp, tauP, Kpva, Kpp, tauL, Kdvdv)
            LHS  = obj.LHSS;
            nVel = size(C, 1);

            tauP = tauP(:);
            tauS = tauS(:);
            tauL = tauL(:);

            Kvva2 = sparse(Kvva2);
            KvS   = Kvva2 * tauS;

            Kavp  = sparse(Kavp);
            KpS   = Kavp * tauS;

            Kpva  = sparse(Kpva);
            KvP   = Kpva * tauP;

            Kpp   = sparse(Kpp);
            KpP   = Kpp * tauP;

            Kdvdv = sparse(Kdvdv);
            KvL   = Kdvdv * tauL;

            LHS(1:nVel, 1:nVel)             = LHS(1:nVel, 1:nVel) + C + KvS + KvL; %+ KvL
            LHS(1:nVel, nVel + 1:end)       = LHS(1:nVel, nVel + 1:end) + KpS;
            LHS(nVel + 1:end, 1:nVel)       = LHS(nVel + 1:end, 1:nVel) + KvP;
            LHS(nVel + 1:end, nVel + 1:end) = LHS(nVel + 1:end, nVel + 1:end) + KpP;
        end

        function tauS1 = computeTauS1(obj,CE,KE)
            normCE = norm(CE,"fro");
            normKE = norm(KE,"fro");
            tauS1  = normCE / normKE;
        end

        function tauS3 = computeTauS3(obj, tauS1, a)
            u   = norm(a,"fro");
            tauS3 = tauS1 ^2 * u^2 / obj.nu;
        end

        function tau = computeTau(obj, tauS1,tauS3)
            tau = (1/(tauS1^obj.residual) + 1/(tauS3^obj.residual)) ^(-1/obj.residual);
        end

        function tauP3 = computeTauP3(obj, tauP1, a, tauS1)
            u   = norm(a,"fro");
            tauP3 = tauP1 * u^2 / obj.nu * tauS1;
        end

    end

end