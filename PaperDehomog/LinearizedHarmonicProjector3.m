classdef LinearizedHarmonicProjector3 < handle

    properties (Access = public)

    end

    properties (Access = private)
        eta
        internalDOFs
        massMatrixBB
        massMatrixGG
        stiffnessMatrix
        fB
        fS
        fG
        perimeter
    end

    properties (Access = private)
        mesh
        boundaryNodes
        density
    end

    methods (Access = public)

        function obj = LinearizedHarmonicProjector3(cParams)
            obj.init(cParams);
            obj.initializeFunctions();
            obj.eta = (10*obj.mesh.computeMeanCellSize)^2;  
            obj.perimeter = obj.density.*(1-obj.density);

            obj.createInternalDOFs();
            obj.computeAllMassMatrix();
            obj.computeStiffnessMatrix();
        end

        function b = solveProblem(obj,bBar,b)
            RHS = obj.computeRHS(bBar);
            LHS = obj.computeLHS(b);
            nInt = size(obj.internalDOFs,2);
            x = [b.fValues(:);zeros(nInt,1);zeros(obj.fG.nDofs,1)];
            res = norm(LHS*x - RHS)/norm(x);
            [resL,resH,resB,resG] = obj.evaluateResidualNorms(bBar,b);
            i = 1;
            theta = 0.5;
            while res(i) > 1e-6
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
            resLnorm = Norm(resL,'L2');
            resHnorm = Norm(resH,'L2');
            resBnorm = Norm(resB,'L2');
            resGnorm = Norm(resG,'L2');
        end

        function b = createVectorFromSolution(obj,x)
            nB = obj.fB.nDofs;
            bV = x(1:2*nB);
            s.fValues = reshape(bV,[],2);
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            b = LagrangianFunction(s);
        end

        function [resL,resH,resB,resG] = evaluateAllResiduals(obj,bBar,b)
            resL = obj.evaluateLossResidual(bBar,b);
            resH = obj.evaluateHarmonicResidual(b);
            resB = obj.evaluateUnitNormResidual(b);
            resG = obj.evaluteGradientNorm(b);
        end

        function difB = evaluateLossResidual(obj,bBar,b)
            difB = DDP(b - bBar,b - bBar);
        end

        function bR = createReshapedFunction(obj,b)
            bR = DomainFunction.create(@(x) reshape(b.evaluate(x),[1 size(b.evaluate(x))]),b.mesh,b.ndimf); 
        end

        function resH = evaluateHarmonicResidual(obj,b)
            bs = b.getVectorFields;            
            f = (-Grad(bs{1}).*bs{2}+Grad(bs{2}).*bs{1});
            s.mesh = obj.mesh;
            s.quadratureOrder = 4;
            s.type = 'ShapeDerivative';
            test = obj.fG;
            rhs  = RHSIntegrator.create(s);
            rhsV = rhs.compute(f,test);
            rhsV(obj.boundaryNodes) = 0;
            Mgg = obj.massMatrixGG; 
            hf = Mgg\rhsV;
            resH = obj.createP1Function((hf));
        end

        function resB = evaluateUnitNormResidual(obj,b)
            f = @(x) obj.computeUnitNormFunction(b,x);
            resB = DomainFunction.create(f,obj.mesh,1); 

           % Mgg = obj.massMatrixGG; 
           % nBf = (Mgg)\nB;            
           % resB = obj.createP1Function(abs(nBf));
        end

        function resG = evaluteGradientNorm(obj,b)
            bS = b.getVectorFields;
            b1  = bS{1};
            b2  = bS{2};
            grad = norm(Grad(b1).*Grad(b1)+Grad(b2).*Grad(b2),2);
            resG = project(grad,'P1D');
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh             = cParams.mesh;
            obj.boundaryNodes    = cParams.boundaryMesh;
            obj.density          = cParams.density;            
        end

        function initializeFunctions(obj)
            obj.fB = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.fS = LagrangianFunction.create(obj.mesh, 1, 'P1');
            obj.fG = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

        function computeAllMassMatrix(obj)
            obj.massMatrixBB = obj.createMassMatrixWithFunction(obj.perimeter);
            obj.massMatrixGG = obj.computeMassMatrix(obj.fG,obj.fG);
        end

        function M = computeMassMatrix(obj,test,trial)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = test;
            s.trial = trial;
            s.quadratureOrder = 2;
            lhs = LHSIntegrator.create(s);
            M = lhs.compute();
        end

        function nB = computeUnitNormFunction(obj,b,x)
            bV = b.evaluate(x);
            b1V  = (bV(1,:,:));
            b2V  = (bV(2,:,:));
            nB   = b1V.^2 + b2V.^2 -1;
        end

        function computeStiffnessMatrix(obj)
            s.test  = obj.fB;
            s.trial = obj.fB;
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            lhs = LHSIntegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end

        function createInternalDOFs(obj)
            bNodes = obj.boundaryNodes;
            iDOFs  = setdiff(1:obj.fS.nDofs,bNodes);
            obj.internalDOFs = iDOFs;
        end

        function Kf = createStiffNessMatrixWithFunction(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fS;
            s.function = f;
            s.mesh     = obj.mesh;
            s.type  = 'StiffnessMatrixWithFunction';
            lhs = LHSIntegrator.create(s);
            Kf = lhs.compute();
        end

        function f = createP1Function(obj,fV)
            s.fValues = fV;
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            f = LagrangianFunction(s);
        end

        function f = createScalarFunctions(obj,b,dir)
            s.fValues = b.fValues(:,dir);
            s.mesh    = obj.mesh;
            s.order   = 'P1';
            f = LagrangianFunction(s);
        end

        function Nf = createAdvectionMatrixWithFunction(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fS;
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 2;
            s.type  = 'AdvectionMatrixWithFunction';
            lhs = LHSIntegrator.create(s);
            Nf = lhs.compute();
        end

        function Mf = createMassMatrixWithFunction(obj,f)
            s.test     = obj.fB;
            s.trial    = obj.fG;
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 2;
            s.type  = 'MassMatrixWithFunction';
            lhs = LHSIntegrator.create(s);
            Mf = lhs.compute();
        end

        function [Mb1,Mb2] = createMassMatrixWithB(obj,b)
            bs = b.getVectorFields;
            b1  = bs{1};
            b2  = bs{2};
            Mb1  = obj.createMassMatrixWithFunction(b1);
            Mb2  = obj.createMassMatrixWithFunction(b2);
        end

        function [Kb1,Kb2,Nb1,Nb2] = computeHarmonicMatrix(obj,b)
            bs = b.getVectorFields;
            b1  = bs{1};
            b2  = bs{2};
            Kb1 = obj.createStiffNessMatrixWithFunction(b1);
            Kb2 = obj.createStiffNessMatrixWithFunction(b2);
            Nb1 = obj.createAdvectionMatrixWithFunction(b1);
            Nb2 = obj.createAdvectionMatrixWithFunction(b2);
            Kb1 = obj.computeReducedAdvectionMatrix(Kb1);
            Kb2 = obj.computeReducedAdvectionMatrix(Kb2);
            Nb1 = obj.computeReducedAdvectionMatrix(Nb1);
            Nb2 = obj.computeReducedAdvectionMatrix(Nb2);
        end

        function LHS = computeLHS(obj,b)
            Mbb  = obj.massMatrixBB;
            K    = obj.stiffnessMatrix;
            nInt = size(obj.internalDOFs,2);
            [Mb1,Mb2] = obj.createMassMatrixWithB(b);
            [Kb1,Kb2,Nb1,Nb2] = obj.computeHarmonicMatrix(b);
            Z  = sparse(obj.fB.nDofs,obj.fB.nDofs);
            Zh = sparse(nInt,nInt);
            Zsg = sparse(nInt,obj.fG.nDofs);
            Zgg = sparse(obj.fG.nDofs,obj.fG.nDofs);
            A  = Mbb + obj.eta*K;
            LHS = [A          ,         Z,(-Kb2+Nb2),Mb1;       ...
                Z          ,         A, (Kb1-Nb1),Mb2;         ...
                (-Kb2+Nb2)',(Kb1-Nb1)',        Zh,Zsg;...
                Mb1',      Mb2',       Zsg',Zgg];
        end

        function RHS = computeRHS(obj,bBar)
            rhsB = obj.computeRHSCostDerivative(bBar);
            rhsH = zeros(size(obj.internalDOFs,2),1);
            rhsU = obj.computeRHSunitNorm();
            RHS  = [rhsB;rhsH;rhsU];
        end

        function rhsB = computeRHSCostDerivative(obj,bBar)
            s.mesh = obj.mesh;
            s.quadType = 3;
            s.type = 'ShapeFunction';
            test = bBar;
            f    = bBar.*obj.perimeter;
            rhs  = RHSIntegrator.create(s);
            rhsB = rhs.compute(f,test);
            rhsB = reshape(rhsB,2,[])';
            rhsB = rhsB(:);
        end

        function rhsV = computeRHSunitNorm(obj)
            s.mesh = obj.mesh;
            s.quadType = 2;
            s.type = 'ShapeFunction';
            f    = ConstantFunction.create(1,obj.mesh);
            rhs  = RHSIntegrator.create(s);
            rhsV = rhs.compute(f,obj.fG);
        end

        function Ared = computeReducedAdvectionMatrix(obj,A)
            iDOFs = obj.internalDOFs;
            Ared = A(:,iDOFs);
        end


    end

end