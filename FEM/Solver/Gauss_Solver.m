classdef Gauss_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS,mesh,bc)
            normVal = Inf;
            tol = 1e-3;
            n = length(LHS);
            x = RHS;
            iter = 0;
            L = tril(LHS);
            U = triu(LHS,1);
            while normVal>tol
                xold=x;
                x=L\(RHS-U*xold);
                if mod(iter,100) == 0
                    %Gauss_Solver.plotSolution(x,mesh,bc,iter)
                    res = LHS*x - RHS;
                    Gauss_Solver.plotRes(res,mesh,bc,iter)
                end
                %normVal=norm(x-xold);
                normVal = norm(LHS*x - RHS);
                iter = iter+1;
                residu(iter) = normVal;
            end
            %save('residuGaussZeros.mat','residu')
        end

        function plotSolution(x,mesh,bc,iter)
            xFull = bc.reducedToFullVector(x);
            s.fValues = reshape(xFull,2,[])';
            s.mesh = mesh;
            s.fValues(:,end+1) = 0;
            s.ndimf = 3;
            xF = P1Function(s);
            %xF.plot();
            xF.print(['sol',num2str(iter)],'Paraview')
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
            xF.print(['Res',num2str(numItr)],'Paraview')
            fclose('all');
        end
    end

end