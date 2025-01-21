classdef LinearizedHarmonicProjector4 < handle

    properties (Access = public)

    end

    properties (Access = private)
        eta
        epsilon
        massMatrixBB
        massMatrixEE
        massMatrixES
        massMatrixGG
        stiffnessMatrix
        divergenceMatrix
        fB
        fE
        fS
        fG
    end

    properties (Access = private)
        mesh
        boundaryNodes
    end

    methods (Access = public)

        function obj = LinearizedHarmonicProjector4(cParams)
            obj.init(cParams);
            obj.initializeFunctions();
            obj.computeAllMassMatrix();
            obj.computeStiffnessMatrix();
            obj.computeDivergenceMatrix();
        end

        function b = solveProblem(obj,bBar,b)
            RHS = obj.computeRHS(bBar);
            LHS = obj.computeLHS(b);
            x = [b.fValues(:);zeros(obj.fE.nDofs,1);zeros(obj.fS.nDofs,1);zeros(obj.fG.nDofs,1)];
            res = norm(LHS*x - RHS)/norm(x);
            [resL,resH,resB,resG] = obj.evaluateResidualNorms(bBar,b);
            i = 1;
            theta = 0.5;
            while res(i) > 1e-12
                xNew   = LHS\RHS;
                x = theta*xNew + (1-theta)*x;
                b   = obj.createVectorFromSolution(x);
                LHS = obj.computeLHS(b);
                i   = i+1;
                res(i) = norm(LHS*x - RHS)/norm(x);
                [resL(i),resH(i),resB(i),resG(i)] = obj.evaluateResidualNorms(bBar,b);
                disp(['iter ',num2str(i),' residual ',num2str(res(i))])
            end
            figure()
            plot(1:i,log([res; resL; resH; resB; resG]))
            legend('LHS*x-RHS','resDistance','resHarmonic','resUnitBall','resGradient')
        end

        function [resLnorm,resHnorm,resBnorm,resGnorm] = evaluateResidualNorms(obj,bBar,b)
            [resL,resH,resB,resG] = obj.evaluateAllResiduals(bBar,b);
            resLnorm = resL.computeL2norm();
            resHnorm = resH.computeL2norm();
            resBnorm = resB.computeL2norm();
            resGnorm = resG.computeL2norm();
        end

        function b = createVectorFromSolution(obj,x)
            bV = x(1:2*obj.fB.nDofs);
            s.fValues = reshape(bV,[],2);
            s.mesh    = obj.mesh;
            s.oder    = 'P1';
            b = LagrangianFunction(s);
        end

        function [resL,resH,resB,resG] = evaluateAllResiduals(obj,bBar,b)
            resL = obj.evaluateLossResidual(bBar,b);
            resH = obj.evaluateHarmonicResidual(b);
            resB = obj.evaluateUnitNormResidual(b);
            resG = obj.evaluteGradientNorm(b);
        end

        function resL = evaluateLossResidual(obj,bBar,b)
            difB = b.fValues - bBar.fValues;
            resL = obj.createP1Function(abs(difB));
        end

        function resH = evaluateHarmonicResidual(obj,b)
            b1  = obj.createScalarFunctions(b,1);
            b2  = obj.createScalarFunctions(b,2);
            Kb1 = obj.createStiffNessMatrixWithFunction(b1);
            Kb2 = obj.createStiffNessMatrixWithFunction(b2);
            Nb1 = obj.createAdvectionMatrixWithFunctionDerivative(b1);
            Nb2 = obj.createAdvectionMatrixWithFunctionDerivative(b2);
            b1V = b1.fValues;
            b2V = b2.fValues;
            resV = (-Kb2'+Nb2')*b1V + (Kb1'-Nb1')*b2V;
            resV(obj.boundaryNodes) = 0;
            resH = obj.createP1Function(abs(resV));
        end

        function resB = evaluateUnitNormResidual(obj,b)
            bNorm = obj.computeUnitNormFunction(b);
            resB = obj.createP1Function(abs(bNorm));
        end

        function resG = evaluteGradientNorm(obj,b)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATICMASS');
            xV = quad.posgp;
            b1  = obj.createScalarFunctions(b,1);
            b2  = obj.createScalarFunctions(b,2);
            db1 = b1.evaluateGradient(xV);
            db2 = b2.evaluateGradient(xV);
            db1V = db1.fValues;
            db2V = db2.fValues;
            grad = db1V.*db1V + db2V.*db2V;
            db   = sum(grad,1);
            s.fValues = db;
            s.mesh    = obj.mesh;
            s.quadrature = xV;
            resG = FGaussDiscontinuousFunction(s);
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.boundaryNodes    = cParams.boundaryMesh;
            obj.epsilon          = cParams.epsilon;
            obj.eta     = (5*obj.mesh.computeMeanCellSize)^2;
        end

        function initializeFunctions(obj)
            obj.fB = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.fE = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            obj.fS = LagrangianFunction.create(obj.mesh,obj.mesh.ndim,'P1');
            obj.fG = LagrangianFunction.create(obj.mesh,1,'P0');
        end

        function computeAllMassMatrix(obj)
            obj.massMatrixBB = obj.computeMassMatrix(obj.fB,obj.fB);
            obj.massMatrixEE = obj.computeMassMatrix(obj.fE,obj.fE);
            obj.massMatrixES = obj.computeMassMatrix(obj.fE,obj.fS);
            obj.massMatrixGG = obj.computeMassMatrix(obj.fG,obj.fG);
        end

        function M = computeMassMatrix(obj,test,trial)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = test;
            s.trial = trial;
            s.quadratureOrder = 'QUADRATICMASS';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();
        end

        function nB = computeUnitNormFunction(obj,b)
            b1V  = b.fValues(:,1);
            b2V  = b.fValues(:,2);
            nB   = b1V.^2 + b2V.^2 -1;
        end

        function computeStiffnessMatrix(obj)
            s.test  = obj.fB;
            s.trial = obj.fB;
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end

        function computeDivergenceMatrix(obj)
            s.test  = obj.fE;
            s.trial = obj.fE;
            s.mesh  = obj.mesh;
            s.type   = 'DivergenceMatrix';
            lhs = LHSintegrator.create(s);
            D = lhs.compute();
            obj.divergenceMatrix = D;
        end

        function Kf = createStiffNessMatrixWithFunction(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fB;
            s.function = f;
            s.mesh     = obj.mesh;
            s.type  = 'StiffnessMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            Kf = lhs.compute();
        end

        function f = createP1Function(obj,fV)
            s.fValues = fV;
            s.order   = 'P1';
            s.mesh    = obj.mesh;
            f = LagrangianFunction(s);
        end

        function f = createScalarFunctions(obj,b,dir)
            s.fValues = b.fValues(:,dir);
            s.order   = 'P1';
            s.mesh    = obj.mesh;
            f = LagrangianFunction(s);
        end

        function Nf = createAdvectionMatrixWithFunctionDerivative(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fB;
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATICMASS';
            s.type  = 'AdvectionMatrixWithFunctionDerivative';
            lhs = LHSintegrator.create(s);
            Nf = lhs.compute();
        end

        function Nf = createAdvectionMatrixWithFunction(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fS;
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATICMASS';
            s.type  = 'AdvectionMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            Nf = lhs.compute();
        end

        function Mf = createMassMatrixWithFunctionDerivative(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fS;
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATICMASS';
            s.type  = 'MassMatrixWithFunctionDerivative';
            lhs = LHSintegrator.create(s);
            Mf = lhs.compute();
        end

        function Mf = createMassMatrixWithFunction(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fG;
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATICMASS';
            s.type  = 'MassMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            Mf = lhs.compute();
        end

        function [Mdb1,Mdb2] = createMassMatrixWithDb(obj,b)
            b1  = obj.createScalarFunctions(b,1);
            b2  = obj.createScalarFunctions(b,2);
            Mdb1  = obj.createMassMatrixWithFunctionDerivative(b1);
            Mdb2  = obj.createMassMatrixWithFunctionDerivative(b2);
        end

        function [Mb1,Mb2] = createMassMatrixWithB(obj,b)
            b1  = obj.createScalarFunctions(b,1);
            b2  = obj.createScalarFunctions(b,2);
            Mb1  = obj.createMassMatrixWithFunction(b1);
            Mb2  = obj.createMassMatrixWithFunction(b2);
        end

        function [Nb1,Nb2] = computeAdvectionWithB(obj,b)
            b1  = obj.createScalarFunctions(b,1);
            b2  = obj.createScalarFunctions(b,2);
            Nb1 = obj.createAdvectionMatrixWithFunction(b1);
            Nb2 = obj.createAdvectionMatrixWithFunction(b2);
        end

        function LHS = computeLHS(obj,b)
            Mbb  = obj.massMatrixBB;
            Mes  = obj.massMatrixES;
            Mee  = 0*obj.massMatrixEE;
            K  = obj.stiffnessMatrix;
            D  = obj.divergenceMatrix;
            De = obj.epsilon*D;
            [Mdb1,Mdb2] = obj.createMassMatrixWithDb(b);
            [Mb1,Mb2] = obj.createMassMatrixWithB(b);
            [Nb1,Nb2] = obj.computeAdvectionWithB(b);

            Zbb = sparse(obj.fB.nDofs,obj.fB.nDofs);
            Zss = sparse(obj.fS.nDofs,obj.fS.nDofs);
            Zgg = sparse(obj.fG.nDofs,obj.fG.nDofs);
            Zeb = sparse(obj.fE.nDofs,obj.fB.nDofs);
            Zsg = sparse(obj.fS.nDofs,obj.fG.nDofs);
            Zeg = sparse(obj.fE.nDofs,obj.fG.nDofs);

            A  = Mbb + obj.eta*K;
            LHS = [          A,       Zbb,        Zeb', (-Nb2+Mdb2), Mb1;...
                Zbb,         A,        Zeb',  (Nb1-Mdb1), Mb2;...
                Zeb,       Zeb,    2*De+2*Mee,       -2*Mes, Zeg;...
                (-Nb2+Mdb2)',(Nb1-Mdb1)',    -2*Mes',         Zss, Zsg;...
                Mb1'        ,       Mb2',       Zeg',        Zsg', Zgg];
        end

        function RHS = computeRHS(obj,bI)
            Mg    = obj.massMatrixGG;
            Mb    = obj.massMatrixBB;
            bI1  = obj.createScalarFunctions(bI,1);
            bI2  = obj.createScalarFunctions(bI,2);
            Ze   = zeros(obj.fE.nDofs,1);
            Zs   = zeros(obj.fS.nDofs,1);
            I    = ones(obj.fG.nDofs,1);
            RHS  = [Mb*bI1.fValues;Mb*bI2.fValues;Ze;Zs;Mg*I];
        end


    end

end