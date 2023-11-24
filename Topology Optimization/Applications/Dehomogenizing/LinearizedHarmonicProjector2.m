classdef LinearizedHarmonicProjector2 < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        eta
        epsilon
        internalDOFs
        massMatrix
        stiffnessMatrix
    end
    
    properties (Access = private)
        mesh
        boundaryNodes
    end
    
    methods (Access = public)
        
        function obj = LinearizedHarmonicProjector2(cParams)
            obj.init(cParams)
            obj.createInternalDOFs();
            obj.computeMassMatrix();
            obj.computeStiffnessMatrix();
        end

        function b = solveProblem(obj,bBar,b)
            RHS = obj.computeRHS(bBar);
            LHS = obj.computeLHS(b);    
            x = [b.fValues(:);zeros(size(obj.internalDOFs,2),1)];
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
            plot(1:i,([res; resL; resH; resB; resG]))
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
            nB = length(x) - length(obj.internalDOFs);
            bV = x(1:nB);
            s.fValues = reshape(bV,[],2);
            s.mesh    = obj.mesh;
            b = P1Function(s);
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
            Nb1 = obj.createAdvectionMatrixWithFunction(b1);
            Nb2 = obj.createAdvectionMatrixWithFunction(b2);
            b1V = b1.fValues;
            b2V = b2.fValues;
            resV = (-Kb2'+Nb2')*b1V + (Kb1'-Nb1')*b2V;
            resV(obj.boundaryNodes) = 0;
            resH = obj.createP1Function(abs(resV));   
% 
%             quad = Quadrature.set(obj.mesh.type);            
%             quad.computeQuadrature('QUADRATICMASS');
% 
%             bG  = b.evaluate(quad.posgp);
%             b1  = bG(1,:,:);
%             b2  = bG(2,:,:);
%             gbG = b.computeGradient(quad);
%             gb1 = gbG.fValues(1:2,:,:);
%             gb2 = gbG.fValues(3:4,:,:);
% 
%             s.fValues(1,:,:) = abs(-b2.*gb1(1,:,:) + b1.*gb2(1,:,:));
%             s.fValues(2,:,:) = abs(-b2.*gb1(2,:,:) + b1.*gb2(2,:,:));
%             s.mesh = obj.mesh;
%             s.quadrature = quad;
%             resH = FGaussDiscontinuousFunction(s);            
        end

        function resB = evaluateUnitNormResidual(obj,b)
            nB   = obj.computeUnitNormFunction(b);
            resB = obj.createP1Function(abs(nB));                         
        end

        function resG = evaluteGradientNorm(obj,b)
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATICMASS');    
            b1  = obj.createScalarFunctions(b,1);
            b2  = obj.createScalarFunctions(b,2);            
            db1 = b1.computeGradient(quad);
            db2 = b2.computeGradient(quad);   
            db1V = db1.fValues;
            db2V = db2.fValues;
            grad = db1V.*db1V + db2V.*db2V;
            db   = sum(grad,1);
            s.fValues = db;
            s.mesh    = obj.mesh;
            s.quadrature = quad;
            resG = FGaussDiscontinuousFunction(s);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh             = cParams.mesh;
           obj.boundaryNodes    = cParams.boundaryMesh;
           obj.epsilon          = cParams.epsilon;
           obj.eta     = (80*obj.mesh.computeMeanCellSize)^2;                            
        end
        
        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();           
            obj.massMatrix = M;
        end    

        function nB = computeUnitNormFunction(obj,b)
            b1V  = b.fValues(:,1);
            b2V  = b.fValues(:,2);
            nB   = b1V.^2 + b2V.^2 -1;
        end

        function computeStiffnessMatrix(obj)        
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.mesh         = obj.mesh;
            s.type         = 'StiffnessMatrix';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end   
        
        function createInternalDOFs(obj)
           bNodes = obj.boundaryNodes;
           iDOFs  = setdiff(1:obj.mesh.nnodes,bNodes); 
           obj.internalDOFs = iDOFs;
        end   

        function Kf = createStiffNessMatrixWithFunction(obj,f)
            s.test     = P1Function.create(obj.mesh, 1);
            s.trial    = P1Function.create(obj.mesh, 1);
            s.function = f;
            s.mesh     = obj.mesh;
            s.type  = 'StiffnessMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            Kf = lhs.compute();
        end

        function f = createP1Function(obj,fV)
           s.fValues = fV;
           s.mesh    = obj.mesh;
           f = P1Function(s);
        end

        function f = createScalarFunctions(obj,b,dir)
            s.fValues = b.fValues(:,dir);
            s.mesh    = obj.mesh;
            f = P1Function(s);
        end

        function Nf = createAdvectionMatrixWithFunction(obj,f)
            s.test     = P1Function.create(obj.mesh, 1);
            s.trial    = P1Function.create(obj.mesh, 1);
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATICMASS';            
            s.type  = 'AdvectionMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            Nf = lhs.compute();
        end        

        function Mf = createMassMatrixWithFunction(obj,f)
            s.test     = P1Function.create(obj.mesh, 1);
            s.trial    = P1Function.create(obj.mesh, 1);
            s.function = f;
            s.mesh     = obj.mesh;
            s.quadratureOrder = 'QUADRATICMASS';            
            s.type  = 'MassMatrixWithFunction';
            lhs = LHSintegrator.create(s);
            Mf = lhs.compute();
        end                

        function Mf = createUnitNormMassMatrix(obj,b)
            nB  = obj.computeUnitNormFunction(b);
            nBF = obj.createP1Function(nB);            
            Mf  = obj.createMassMatrixWithFunction(nBF);
        end

        function [Kb1,Kb2,Nb1,Nb2] = computeHarmonicMatrix(obj,b)
            b1  = obj.createScalarFunctions(b,1);
            b2  = obj.createScalarFunctions(b,2);            
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
            M  = obj.massMatrix;
            K  = obj.stiffnessMatrix;
            Mf = obj.createUnitNormMassMatrix(b);
            [Kb1,Kb2,Nb1,Nb2] = obj.computeHarmonicMatrix(b);
            Z  = sparse(size(b.fValues,1),size(b.fValues,1));
            Zh = sparse(size(obj.internalDOFs,2),size(obj.internalDOFs,2));
            A  = M + obj.eta*K + obj.epsilon*Mf;
            LHS = [A, Z, (-Kb2+Nb2); Z, A, (Kb1-Nb1); (-Kb2+Nb2)' (Kb1-Nb1)' Zh];
        end

        function RHS = computeRHS(obj,bI)
            M    = obj.massMatrix;
            bI1  = obj.createScalarFunctions(bI,1);
            bI2  = obj.createScalarFunctions(bI,2);               
            Z    = zeros(size(obj.internalDOFs,2),1);            
            RHS  = [M*bI1.fValues;M*bI2.fValues;Z];
        end

       function Ared = computeReducedAdvectionMatrix(obj,A)
           iDOFs = obj.internalDOFs;
           Ared = A(:,iDOFs);
       end        
        
        
    end
    
end