classdef conjugateGradient_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS,mesh,bc)
            tol = 1e-6;
            n = length(RHS);
            x = RHS; 
            r = RHS - LHS * x; 
            p = r; 
            rsold = r' * r;
            iter = 0;

            hasNotConverged = true;

            while hasNotConverged
                Ap = LHS * p;
                alpha = rsold / (p' * Ap);
                x = x + alpha * p;
                r = r - alpha * Ap;
                rsnew = r' * r;

                %hasNotConverged = sqrt(rsnew) > tol;
                hasNotConverged = max(LHS*x - RHS) > tol;

                p = r + (rsnew / rsold) * p;
                rsold = rsnew;
                iter = iter + 1;
                residu(iter) = norm(LHS*x - RHS); %Ax - b
                res = LHS*x - RHS;
                
                conjugateGradient_Solver.plotSolution(x,mesh,bc,iter)
                conjugateGradient_Solver.plotRes(res,mesh,bc,iter)
            end
            %save('residuConjugateZeros.mat', 'residu')
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