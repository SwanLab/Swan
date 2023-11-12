classdef preconditionedConjugateGradient < Solver
    
    methods (Static)
        
        function x = solve(LHS,RHS,mesh,bc)
            tol = 1e-6;
            n = length(RHS);
            x = zeros(n,1);
            r = RHS - LHS * x;
            %             z = ModalTesting.applyPreconditioner(M,r);
            z = preconditionedConjugateGradient.applyPreconditioner(r);
            %             z=r-z;
            p = z;
            rzold = r' * z;
            iter = 0;

            hasNotConverged = true;

            while hasNotConverged
                Ap = LHS * p;
                alpha = rzold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                %                 z = ModalTesting.applyPreconditioner(M,r);
                z = preconditionedConjugateGradient.applyPreconditioner(r);
                rznew = r' * z;

                %hasNotConverged = sqrt(rsnew) > tol;
                hasNotConverged = norm(r) > tol;

                p = z + (rznew / rzold) * p;
                rzold = rznew;
                iter = iter + 1;
                residu(iter) = norm(LHS*x - RHS); %Ax - b
                res = LHS*x - RHS;
                
                preconditionedConjugateGradient.plotSolution(x,mesh,bc,iter)
                preconditionedConjugateGradient.plotRes(res,mesh,bc,iter)
            end
            save('residuConjugateZeros.mat', 'residu')
        end
        %
        function z = applyPreconditioner(obj,r)
            lhs=obj.Kmodal;
            phi=obj.eigenVec;
            r1=phi'*r;
            zP=lhs\r1;
            z=phi*zP;
            %z = r;
            z = r-z;
        end
        %
        
        function [basis,basisVec,eigenVec] = computeBasis(obj,K)
            [eigenVec,D]=eigs(K,obj.nbasis,'smallestabs');
            psi = K*eigenVec;

            for i = 1:size(eigenVec,2)
                b = eigenVec(:,i);
                b1 = obj.boundaryConditions.reducedToFullVector(b);
                basis{i} = reshape(b1,2,[])';
                basisVec{i}= b1;
                a = psi(:,i);
                a1=obj.boundaryConditions.reducedToFullVector(a);
                psiD{i}=reshape(a1,2,[])';
                %  bC{i} = b;

                % sF.fValues = reshape(b1,2,[])';
                % sF.mesh    = mesh;
                % bF{i} =P1Function(sF);
                % bF.plot
            end

        end
        
        function plotSolution(x,mesh,bc,numItr)
            xFull = bc.reducedToFullVector(x);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['sol',num2str(numItr)],'GiD')
            fclose('all');
        end
        
        function plotRes(res,mesh,bc,numItr)
            xFull = bc.reducedToFullVector(res);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['Res',num2str(numItr)],'GiD')
            fclose('all');
        end
    end
end

