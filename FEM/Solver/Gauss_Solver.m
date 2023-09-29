classdef Gauss_Solver < Solver

    methods (Static)

        function x = solve(LHS,RHS,mesh,bc)
            normVal = Inf;
            tol = 1e-6;
            n = length(LHS);
            x = rand(n,1);
            numItr = 0;
            L = tril(LHS);
            U = triu(LHS,1);
            while normVal>tol
                xold=x;
                x=L\(RHS-U*xold);
                %Gauss_Solver.plotSolution(x,mesh,bc,numItr)
                normVal=norm(x-xold);
                numItr = numItr+1;
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
            xF.print(['sol',num2str(numItr)],'Paraview')
            fclose('all');
        end
    end

end